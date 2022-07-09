Here are examples of using simulatorBMP() with a simplified 2D model of the vocal tract.
The results are interesting, but not exact. Note that the simulation runs in 2D. A point source in 2D behaves like a cylindrical one in 3D.
The sections shown in "TubeModel_a.bmp" are not cylinder sections. All sections in the z dimension would be exactly the same.
In a cylindrical section, if you double the diameter of a tube, the area increases proportionally to the diameter^2.
In 2D models, the equivalent area grows proportional to the diameter.

To run the model, use simulatorModel.m
The BMP() simulation function directory must be in the Matlab path, unless you put all the files in the same directory.
Select the image to run by changing the imageName variable in the script.
There are two boolean variables in the script.
Use writeAudioFile=true; to create an audio file with the sound results. The file name will be the same as the bmp, but with ".wav"
Use flagGlottal=true to generate smooth waveform glottal pulses (see the quote in the glottalSource.m function).
Use flagGlottal=false to generate a rectangular pulse train.

Jorge Petrosino, june 2022
Universidad Nacional de Lanus
