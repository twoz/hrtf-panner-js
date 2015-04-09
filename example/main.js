window.onload = function () {
	// initialize Audio Context
	try {
		var audioContext = new (window.AudioContext || window.webkitAudioContext)();		
	}
	catch (e) {
		alert("Web Audio API is not supported in this browser");
	}

	// load hrir to the container
	var hrtfContainer = new HRTFContainer();
	hrtfContainer.loadHrir("../hrir/kemar_L.bin");

	// create audio source node from the <audio> element
	var sourceNode = audioContext.createMediaElementSource(document.getElementById("player"));
	var gain = audioContext.createGain();
	gain.gain.value = 0.3;
	sourceNode.connect(gain);

	// create new hrtf panner, source node gets connected automatically
	var panner = new HRTFPanner(audioContext, gain, hrtfContainer);

	// connect the panner to the destination node
	panner.connect(audioContext.destination);

	// animate source
	var t = 0;
	var x, y, z;
	setInterval(function () {
		x = Math.sin(t);
		y = Math.cos(t);
		z = 0;
		t += 0.05;
		var cords = cartesianToInteraural(x, y, z);
		panner.update(cords.azm, cords.elv);
	}, 50);
}