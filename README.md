# ELEC-E5660 Project: LPC-based Vocal Synthesizer
## Student: Francesco Colotti
## Student number: 101613215

## Description
This project aimed at creating a vocal synthesizer based on Linear Predictive
Coding that is able to reproduce basic vowel sounds and allows the user to change some parameters, like an ADSR envelope, in real time. It was coded entirely in MATLAB, with the aid of some external programs to record and clean audio samples, and can be played using an external MIDI keyboard.

The script works with any audio device, though an ASIO device is preferable to obtain a lower latency when playing.

During development and testing the following setup was used:
* Windows 10
* MIDI device: Native Instruments Komplete Kontrol M32
* Audio device: Focusrite 2i2 3rd generation
* Driver: Focusrite USB ASIO Driver (when testing with said audio card) OR FL Studio ASIO (to achieve decent latency even without the USB audio card)

## Setup
To run the code it is necessary to setup it properly in `main.m`. In particular it is vital to correctly write down the MIDI device name (if one wants to use it) which must be the same that gets printed when calling `mididevinfo`. The audio device can be setup explicitely, otherwise the script will use the current default one; the code used in testing has been left in the script to reduce the time needed to select another device (e.g. with an ASIO driver).

## MIDI CCs
* 14: LPC order (10 to 100)
* 15: reverb dry/wet mix
* 18: envelope attack time
* 19: envelope decay time
* 20: envelope sustain level
* 21: envelope release time

If necessary the MIDI CC numbers can be modified inside the synth function.

## Requirements
The following MATLAB toolboxes are required
(version 23.2 was used for everything when developing):
```
Signal Processing Toolbox
DSP System Toolbox
Audio Toolbox
```
The list was obtained by running:
```
[fList,pList] = matlab.codetools.requiredFilesAndProducts('main.m');
```
## Included files
* `main.m`: the main script containing the synthesizer function
* `LPC_coeff_creator.m`: the script used to create the LPC coefficient tables. Must be run before `main.m` if the tables have not been created already
* `report_graphs.m`: the script used to create the graphs in the report
## Folders
* `coeffs`: contains LPC coefficient tables
* `coeffs_N`: contains tables with only LPC coefficients of order N
* `docs`: documentation (report) and content used in the report (e.g. audio files)
* `samples`: samples used to create the LPC coefficients