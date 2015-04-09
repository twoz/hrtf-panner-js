function HRTFContainer() {

	var hrir =  {};
	var triangulation = {
		points: [],
		triangles: []
	};

	this.loadHrir = function(file, onLoad) {
		var oReq = new XMLHttpRequest();
		oReq.open("GET", file, true);
		oReq.responseType = "arraybuffer";
		oReq.onload = function(oEvent) {
			var arrayBuffer = oReq.response;
			if (arrayBuffer) {
				var rawData = new Float32Array(arrayBuffer);
				var ir = {};
				ir.L = {};
				ir.R = {};
				var azimuths = [-90, -80, -65, -55, -45, -40, -35, -30, -25, -20,
					-15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 55, 65, 80, 90];
				var points = [];

				var hrirLength = 200;
				var k = 0;
				for (var i = 0; i < azimuths.length; ++i) {
					azi = azimuths[i];
					ir['L'][azi] = {};
					ir['R'][azi] = {};

					// -90 deg elevation
					ir['L'][azi][-90] = rawData.subarray(k, k + hrirLength);
					k += hrirLength;
					ir['R'][azi][-90] = rawData.subarray(k, k + hrirLength);
					k += hrirLength;

					points.push([azi, -90]);
					// 50 elevations: -45 + 5.625 * (0:49)
					for (var j = 0; j < 50; ++j) {
						var elv = -45 + 5.625 * j;
						ir['L'][azi][elv] = rawData.subarray(k, k + hrirLength);
						k += hrirLength;
						ir['R'][azi][elv] = rawData.subarray(k, k + hrirLength);
						k += hrirLength;
						points.push([azi, elv]);
					}

					// 270 deg elevation
					ir['L'][azi][270] = rawData.subarray(k, k + hrirLength);
					k += hrirLength;
					ir['R'][azi][270] = rawData.subarray(k, k + hrirLength);
					k += hrirLength;
					points.push([azi, 270]);
				}

				hrir = ir;
				triangulation.triangles = Delaunay.triangulate(points);
				triangulation.points = points;
				if (typeof onLoad !== "undefined")
					onLoad();
			}
			else {
				throw new Error("Failed to load HRIR");
			}
		};
		oReq.send(null);
	}

	this.interpolateHRIR = function(azm, elv) {
		var triangles = triangulation.triangles;
		var points = triangulation.points;
		var i = triangles.length - 1;
		var A, B, C, X, T, invT, det, g1, g2, g3;
		while (true) {
			A = points[triangles[i]]; i--;
			B = points[triangles[i]]; i--;
			C = points[triangles[i]]; i--;
			T = [A[0] - C[0], A[1] - C[1],
				 B[0] - C[0], B[1] - C[1]];
			invT = [T[3], -T[1], -T[2], T[0]];
			det = 1 / (T[0] * T[3] - T[1] * T[2]);
			for (var j = 0; j < invT.length; ++j)
				invT[j] *= det;
			X = [azm - C[0], elv - C[1]];
			g1 = invT[0] * X[0] + invT[2] * X[1];
			g2 = invT[1] * X[0] + invT[3] * X[1];
			g3 = 1 - g1 - g2;
			if (g1 >= 0 && g2 >= 0 && g3 >= 0) {
				var hrirL = new Float32Array(200);
				var hrirR = new Float32Array(200);
				for (var i = 0; i < 200; ++i) {
					hrirL[i] = g1 * hrir['L'][A[0]][A[1]][i] +
						g2 * hrir['L'][B[0]][B[1]][i] +
						g3 * hrir['L'][C[0]][C[1]][i];
					hrirR[i] = g1 * hrir['R'][A[0]][A[1]][i] +
						g2 * hrir['R'][B[0]][B[1]][i] +
						g3 * hrir['R'][C[0]][C[1]][i];
				}
				return [hrirL, hrirR];
			}
			else if (i < 0) {
				break;
			}
		}
		return [new Float32Array(200), new Float32Array(200)];
	}
}

/*
	audioContext - an instance of the Web Audio API AudioContext
	sourceNode - any instance of the AudioNode
	hrtfContainer - an instance of HRTFContainer class this panner will use
	Note that for this panner to work properly the source must be MONO.
*/
function HRTFPanner(audioContext, sourceNode, hrtfContainer) {

	function HRTFConvolver() {
		this.buffer = audioContext.createBuffer(2, 200, audioContext.sampleRate);
		this.convolver = audioContext.createConvolver();
		this.convolver.normalize = false;
		this.convolver.buffer = this.buffer;
		this.gainNode = audioContext.createGain();

		this.convolver.connect(this.gainNode);

		this.fillBuffer = function(hrirLR) {
			var bufferL = this.buffer.getChannelData(0);
			var bufferR = this.buffer.getChannelData(1);
			for (var i = 0; i < this.buffer.length; ++i) {
				bufferL[i] = hrirLR[0][i];
				bufferR[i] = hrirLR[1][i];
			}
			this.convolver.buffer = this.buffer;
		}
	}

	var currentConvolver = new HRTFConvolver();
	var targetConvolver = new HRTFConvolver();

	var loPass = audioContext.createBiquadFilter();
	var hiPass = audioContext.createBiquadFilter();
	loPass.type = "lowpass";
	loPass.frequency.value = 200;
	hiPass.type = "highpass";
	hiPass.frequency.value = 200;

	var source = sourceNode;
	source.channelCount = 1;
	source.connect(loPass);
	source.connect(hiPass);
	hiPass.connect(currentConvolver.convolver);
	hiPass.connect(targetConvolver.convolver);


	/* 
		Connects this panner to the destination node.
	*/
	this.connect = function(destination) {
		loPass.connect(destination);
		currentConvolver.gainNode.connect(destination);
		targetConvolver.gainNode.connect(destination);
	}

	/* 
		Connects a new source to this panner and disconnects the previous one.
	*/
	this.setSource = function(newSource) {
		source.disconnect(loPass);
		source.disconnect(hiPass);
		newSource.connect(loPass);
		newSource.connect(hiPass);
		source = newSource;
	}

	/*
		Sets a cut-off frequency below which input signal won't be spatialized.
	*/
	this.setCrossoverFrequency = function(freq) {
		loPass.frequency.value = freq;
		hiPass.frequency.value = freq;
	}


	/*
		Updates the current Head Related Impulse Response.
		Azimuth and elevation are coordinates of the source in the Interaural-Polar 
		coordinate system RELATIVE to the listener.
		This is supposed to be called each time a listener or source position changes.	
	*/
	this.update = function(azimuth, elevation) {
		targetConvolver.fillBuffer(hrtfContainer.interpolateHRIR(azimuth, elevation));
		// start crossfading
		var crossfadeDuration = 25;
		targetConvolver.gainNode.gain.setValueAtTime(0, audioContext.currentTime);
		targetConvolver.gainNode.gain.linearRampToValueAtTime(1,
			audioContext.currentTime + crossfadeDuration / 1000);
		currentConvolver.gainNode.gain.setValueAtTime(1, audioContext.currentTime);
		currentConvolver.gainNode.gain.linearRampToValueAtTime(0,
			audioContext.currentTime + crossfadeDuration / 1000);
		// swap convolvers
		var t = targetConvolver;
		targetConvolver = currentConvolver;
		currentConvolver = t;
	}

}







// ------------------------   HELPER FUNCTIONS  -------------------------------


/* 
	Params:
	x1 - axis that passes through the ears from left to right
	x2 - axis that passes "between the eyes" and points ahead
	x3 - axis that points "up"

	Returns point in interaural coordinates.
*/
function cartesianToInteraural(x1, x2, x3) {
	var r = Math.sqrt(x1 * x1 + x2 * x2 + x3 * x3);
	var azm = rad2deg(Math.asin(x1 / r));
	var elv = rad2deg(Math.atan2(x3, x2));
	if (x2 < 0 && x3 < 0)
		elv += 360;
	return { r: r, azm: azm, elv: elv };
}

function interauralToCartesian(r, azm, elv) {
	azm = deg2rad(azm);
	elv = deg2rad(elv);
	var x1 = r * Math.sin(azm);
	var x2 = r * Math.cos(azm) * Math.cos(elv);
	var x3 = r * Math.cos(azm) * Math.sin(elv);
	return { x1: x1, x2: x2, x3: x3 };
}

function deg2rad(deg) {
	return deg * Math.PI / 180;
}

function rad2deg(rad) {
	return rad * 180 / Math.PI;
}
// ----------------------------------------------------------------------------




// ---------------------------   DELAUNAY   -----------------------------------
// Delaunay triangulation for the hrir interpolation algorithm.
// Delaunay.js - code by ironwallaby (https://github.com/ironwallaby/delaunay)

var Delaunay;

(function() {
	"use strict";

	var EPSILON = 1.0 / 1048576.0;

	function supertriangle(vertices) {
		var xmin = Number.POSITIVE_INFINITY,
			ymin = Number.POSITIVE_INFINITY,
			xmax = Number.NEGATIVE_INFINITY,
			ymax = Number.NEGATIVE_INFINITY,
			i, dx, dy, dmax, xmid, ymid;

		for (i = vertices.length; i--;) {
			if (vertices[i][0] < xmin) xmin = vertices[i][0];
			if (vertices[i][0] > xmax) xmax = vertices[i][0];
			if (vertices[i][1] < ymin) ymin = vertices[i][1];
			if (vertices[i][1] > ymax) ymax = vertices[i][1];
		}

		dx = xmax - xmin;
		dy = ymax - ymin;
		dmax = Math.max(dx, dy);
		xmid = xmin + dx * 0.5;
		ymid = ymin + dy * 0.5;

		return [
		  [xmid - 20 * dmax, ymid - dmax],
		  [xmid, ymid + 20 * dmax],
		  [xmid + 20 * dmax, ymid - dmax]
		];
	}

	function circumcircle(vertices, i, j, k) {
		var x1 = vertices[i][0],
			y1 = vertices[i][1],
			x2 = vertices[j][0],
			y2 = vertices[j][1],
			x3 = vertices[k][0],
			y3 = vertices[k][1],
			fabsy1y2 = Math.abs(y1 - y2),
			fabsy2y3 = Math.abs(y2 - y3),
			xc, yc, m1, m2, mx1, mx2, my1, my2, dx, dy;

		/* Check for coincident points */
		if (fabsy1y2 < EPSILON && fabsy2y3 < EPSILON)
			throw new Error("Eek! Coincident points!");

		if (fabsy1y2 < EPSILON) {
			m2 = -((x3 - x2) / (y3 - y2));
			mx2 = (x2 + x3) / 2.0;
			my2 = (y2 + y3) / 2.0;
			xc = (x2 + x1) / 2.0;
			yc = m2 * (xc - mx2) + my2;
		}

		else if (fabsy2y3 < EPSILON) {
			m1 = -((x2 - x1) / (y2 - y1));
			mx1 = (x1 + x2) / 2.0;
			my1 = (y1 + y2) / 2.0;
			xc = (x3 + x2) / 2.0;
			yc = m1 * (xc - mx1) + my1;
		}

		else {
			m1 = -((x2 - x1) / (y2 - y1));
			m2 = -((x3 - x2) / (y3 - y2));
			mx1 = (x1 + x2) / 2.0;
			mx2 = (x2 + x3) / 2.0;
			my1 = (y1 + y2) / 2.0;
			my2 = (y2 + y3) / 2.0;
			xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
			yc = (fabsy1y2 > fabsy2y3) ?
			  m1 * (xc - mx1) + my1 :
			  m2 * (xc - mx2) + my2;
		}

		dx = x2 - xc;
		dy = y2 - yc;
		return { i: i, j: j, k: k, x: xc, y: yc, r: dx * dx + dy * dy };
	}

	function dedup(edges) {
		var i, j, a, b, m, n;

		for (j = edges.length; j;) {
			b = edges[--j];
			a = edges[--j];

			for (i = j; i;) {
				n = edges[--i];
				m = edges[--i];

				if ((a === m && b === n) || (a === n && b === m)) {
					edges.splice(j, 2);
					edges.splice(i, 2);
					break;
				}
			}
		}
	}

	Delaunay = {
		triangulate: function(vertices, key) {
			var n = vertices.length,
				i, j, indices, st, open, closed, edges, dx, dy, a, b, c;

			/* Bail if there aren't enough vertices to form any triangles. */
			if (n < 3)
				return [];

			/* Slice out the actual vertices from the passed objects. (Duplicate the
			 * array even if we don't, though, since we need to make a supertriangle
			 * later on!) */
			vertices = vertices.slice(0);

			if (key)
				for (i = n; i--;)
					vertices[i] = vertices[i][key];

			/* Make an array of indices into the vertex array, sorted by the
			 * vertices' x-position. */
			indices = new Array(n);

			for (i = n; i--;)
				indices[i] = i;

			indices.sort(function(i, j) {
				return vertices[j][0] - vertices[i][0];
			});

			/* Next, find the vertices of the supertriangle (which contains all other
			 * triangles), and append them onto the end of a (copy of) the vertex
			 * array. */
			st = supertriangle(vertices);
			vertices.push(st[0], st[1], st[2]);

			/* Initialize the open list (containing the supertriangle and nothing
			 * else) and the closed list (which is empty since we havn't processed
			 * any triangles yet). */
			open = [circumcircle(vertices, n + 0, n + 1, n + 2)];
			closed = [];
			edges = [];

			/* Incrementally add each vertex to the mesh. */
			for (i = indices.length; i--; edges.length = 0) {
				c = indices[i];

				/* For each open triangle, check to see if the current point is
				 * inside it's circumcircle. If it is, remove the triangle and add
				 * it's edges to an edge list. */
				for (j = open.length; j--;) {
					/* If this point is to the right of this triangle's circumcircle,
					 * then this triangle should never get checked again. Remove it
					 * from the open list, add it to the closed list, and skip. */
					dx = vertices[c][0] - open[j].x;
					if (dx > 0.0 && dx * dx > open[j].r) {
						closed.push(open[j]);
						open.splice(j, 1);
						continue;
					}

					/* If we're outside the circumcircle, skip this triangle. */
					dy = vertices[c][1] - open[j].y;
					if (dx * dx + dy * dy - open[j].r > EPSILON)
						continue;

					/* Remove the triangle and add it's edges to the edge list. */
					edges.push(
					  open[j].i, open[j].j,
					  open[j].j, open[j].k,
					  open[j].k, open[j].i
					);
					open.splice(j, 1);
				}

				/* Remove any doubled edges. */
				dedup(edges);

				/* Add a new triangle for each edge. */
				for (j = edges.length; j;) {
					b = edges[--j];
					a = edges[--j];
					open.push(circumcircle(vertices, a, b, c));
				}
			}

			/* Copy any remaining open triangles to the closed list, and then
			 * remove any triangles that share a vertex with the supertriangle,
			 * building a list of triplets that represent triangles. */
			for (i = open.length; i--;)
				closed.push(open[i]);
			open.length = 0;

			for (i = closed.length; i--;)
				if (closed[i].i < n && closed[i].j < n && closed[i].k < n)
					open.push(closed[i].i, closed[i].j, closed[i].k);

			/* Yay, we're done! */
			return open;
		},
		contains: function(tri, p) {
			/* Bounding box test first, for quick rejections. */
			if ((p[0] < tri[0][0] && p[0] < tri[1][0] && p[0] < tri[2][0]) ||
			   (p[0] > tri[0][0] && p[0] > tri[1][0] && p[0] > tri[2][0]) ||
			   (p[1] < tri[0][1] && p[1] < tri[1][1] && p[1] < tri[2][1]) ||
			   (p[1] > tri[0][1] && p[1] > tri[1][1] && p[1] > tri[2][1]))
				return null;

			var a = tri[1][0] - tri[0][0],
				b = tri[2][0] - tri[0][0],
				c = tri[1][1] - tri[0][1],
				d = tri[2][1] - tri[0][1],
				i = a * d - b * c;

			/* Degenerate tri. */
			if (i === 0.0)
				return null;

			var u = (d * (p[0] - tri[0][0]) - b * (p[1] - tri[0][1])) / i,
				v = (a * (p[1] - tri[0][1]) - c * (p[0] - tri[0][0])) / i;

			/* If we're outside the tri, fail. */
			if (u < 0.0 || v < 0.0 || (u + v) > 1.0)
				return null;

			return [u, v];
		}
	};

	if (typeof module !== "undefined")
		module.exports = Delaunay;
})();