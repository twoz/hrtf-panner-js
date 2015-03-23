USAGE:

/* This will load hrirs */
var hrtfContainer = new HRTFContainer("path_to_hrir_folder");

var panner = new HRTFPanner(audioContext, hrtfContainer, sourceNode);

panner.update(azimuth, elevation);