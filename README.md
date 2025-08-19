# UNDT-Hydrophone

This repo contains info on how to use the UNDT Lab's Hydrophone rig.

## Hardware Instructions

### User Manual

Our system is the Precision Acoustic Fibre-Optic Hydrophone Version 1 (See manual in _Doccumentation_, however the Version 2 model is VERY similar and its manual better explains the cleaning procedure: [Fibre-Optic Hydrophone V2 User Guide](https://www.acoustics.co.uk/wp-content/uploads/2024/02/Fibre-optic-Hydrophone-System-Version-2-User-Guide-v2.0-1.pdf). A paper copy of the V1 instructions should be stored with the hardware.

 Any notable differences to the V2 system are explained below:

- The optical fibre connector is different. Our version must be oriented correctly then screwed in. It is also cleaned differently:
	- To clean the optical end face wipe it along the cloth of the Cletop-S fibre connector cleaner. The end face has and angle to it so you must  find the correct orientation by feel.
	- It can be cleaned with the Neoclean but the Cletop is more convenient.
- Our hardware runs on an older software. See the README for advice on how to install this.
	- FOH Control must be run with admin privilages.	
	- If all else fails and the FOH Control won’t install it is possible to update the software to the current version (for a price that was too steep at the time of writing).
- Our older connector uses different adaptor tips on the EasyGet2. Use ….. to inspect the fibre, and ….  to inspect the socket.

### Bonus Tips
- Use clean water: De-Ionised or Distilled. This reduces the risk of getting particulates/calcification on the sensor.
- Make sure no bubbles are stuck to the sensor tip.
- Check the ITF to ensure that the sensor is functioning correctly.
- If in doubt email tech support they have been quite helpful!
- Have patience when coiling up the optical fibre. Roll it between your fingers so it doesn’t get twisted.
- Minimise the number of times you need to clean the optical connections by leaving the fibre connected between experiments. However, make sure the sensor is protected whilst it is out.
- Instead of pulling the sensor cover off the fibre optic, slide it backwards. This reduces the risk of damaging the sensor tip when reinstalling the cover.

### Equipment List
- Hydrophone (FOH)
	- USB cable
	- Power cable
- HandyScope
- BNC cables x2
- Zaber stage rig
	- USB 
	- Power cable
- FOH Sensors
	- Sensor mount
- Dimension EasyGet endface inspector
- Neoclean-E fiber connector cleaner
- Cletop-S fiber connector cleaner

Please use the [Hardware Log](https://uob-my.sharepoint.com/:x:/r/personal/gv19838_bristol_ac_uk/Documents/PhD/Hydrophone/UNDT-Hydrophone/Hardware%20Log.xlsx?d=wee4d49348d0c4a80830cae2a6db38e84&csf=1&web=1&e=UEH9Df) to note any changes or damage to the apparatus to avoid people accidentally making uncalibrated measurements.

### Dependencies
_Zaber Launcher_ (or _Zaber Console_ if having compatibility problems): https://www.zaber.com/software

_Zaber COM Driver_ (If USB port not recognised check that you have the appropriate FTDI USB driver installed for your connection.): https://www.ftdichip.com/Support/Documents/InstallGuides.htm 

_Zaber Motion Control_: https://software.zaber.com/motion-library/docs/tutorials/install/matlab

_TiePie Drivers_: [https://www.tiepie.com/en/download/archive](https://www.tiepie.com/en/download/archive). "TiePie USB driver 8.1.9"

LibTiePie may need this Matlab AddOn: MATLAB Support for MinGW-w64 C/C++/Fortran Compiler

_FOH Control_: Installed with wizard over USB from the hardware. 

This can be quite difficult to get running...

Does not work without the Microsoft .net Framework Version 3.5 Runtime.
https://dotnet.microsoft.com/en-us/download/dotnet-framework

FOH Control must be run with Admin privilages.


Matlab (2022a) Addons:

_MATLAB Support for MinGW-w64 C/C++/Fortran Compiler_ 

_Matlab LibTiePie 0.6.4 instrument driver for USB scopes_

## Scan and Analysis Programs
Scanning and analysis programs are seperated so that the raw data is preserved and analysis methods can be independantly developed without having to repeat scans. However, there are still different types of scan e.g. ome scans use a different rasterisation method. This means that scan and analysis programs must be matched. 

All programs save the data and scan paremeters to the DataOut folder.

### Scan Programs:

**Standard Scans:**

Note that these are underdeveloped as they are rarely used at the moment.

_3D:_ Cuboid scan volume. Arbitrary voxel resolution. Can be slow.

**Triggered Scans:**

These progams wait for a trigger input at each point in the raster before taking a measurement.

_LiveWaveformView:_ No rasterisation. Just a live view of the signal measured by the FOH.

_3D_WithTrigger:_ Cuboid scan volume. Arbitrary voxel resolution. Can be slow.

_3D_Parameterised:_ All scan parameters have been moved to the _Input_ file to make it easier to make changes.

_TriPlanar:_ Three orthogonal scan planes of independant size. This is a quick way of inspecting a large volume with relatively high resolution.

**Auxiliary Scans:**

_ScanVolumeChecker:_ Quick program to inspect and test run the scan volume. No data is recorded.

### Analysis Programs

_Blank:_ Loads data in. Get creative!

_TriPlanar:_ Plots peak voltages in 3D.

_1/2/3D:_ Plots data from equivalent scans. A little bit inconsistent. Best used as an example.

_LiveWaveformAnalysis:_ Plots some confusing stuff, but also a useful plot of the waveform with a slider to scroll through the time axis. Could do with improving.

