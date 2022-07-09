# simulateBMP v3.0
MATLAB Function based on k-Wave to simulate 2D acoustic wave propagation in air (audible range)
The acoustic space is defined by a picture (BMP file) using color codes to define sensors, sources, reflectant surfaces, and absorption materials. It was designed for maximum simplicity, but lets users modify multiple parameters such as the source signal, including múltiple sources with different signal definition, and the use of an audio signal to drive the sources.
It is an improved version of simulateImage256() presented in ICA2019 (International Congress on Acoustics) on Aachen, Germany.
"MATLAB-based simulation software as teaching aid for physical acoustics"
Jorge Petrosino, Lucas Landini, Georgina Lizaso, Ian Kuri, Ianina Canalis
Universidad Nacional de Lanus, Argentina
23rd International Congress on Acoustics, 2019.
http://pub.dega-akustik.de/ICA2019/data/articles/001342.pdf
This youTube video shows how it can be used
https://youtu.be/aKQ2G7RMh64

Designed to work with the k-Wave toolbox, available for free at: 
http://www.k-wave.org/

Installing kwave is very simple. After download the toolbox, include its directory to de Matlab pass (with subdirectories).

--------------------------
To run samples, use helloWord.m

--------------------------
There are also samples using vocal-tract models in 2D. It runs the simulation, shows graphs of the signal and its spectrum and plays the sound.
Use simulateModel.m with imageName='TubeModel_a.bmp', for example, and select writeAudioFile=true, for save audio file. Change imageName to 'TubeModel_e.bmp' or 'TubeModel_o.bmp' to compare results.

Jorge Petrosino, june 2022
Universidad Nacional de Lanus, Argentina
