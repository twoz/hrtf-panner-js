## HRTF Panner

#####Binaural panner using Head-Related Transfer Functions to convert mono audio source to 3D stereo.

###USAGE:


```javascript
// load Head-Related-Impulse-Response from the file
hrtfContainer.loadHrir("path_to_hrir_file");
// create a new panner, source is automatically connected
var panner = new HRTFPanner(audioContext, sourceNode);
// connect the panner to the destination node
panner.connect(audioContext.destination);

// call this each time the relative position between the observer and listener changes
panner.update(azimuth, elevation);
```

