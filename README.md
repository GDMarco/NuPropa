<p align="right">
  <img src="images/nupropaLogo.png"
       alt="nupropa logo"
       width="300">
</p>

# `νpropa`

Plugin to the CRPropa code to implement the propagation of highly energetic neutrinos through photon and neutrino backgrounds, e.g. the cosmic microwave and neutrino backgrounds. The production of the secondaries of interaction is treated.

(data/ folder contains the interaction rates employed in the interaction modules. Its position has to be sync in the code.)

<!-- ## Scientific scenarios -->

## Installation 

To install `νpropa`, you will need to have CRPropa 3 installed (check [CRPropa](https://github.com/CRPropa/CRPropa3) for the latest version).

Steps to install this plugin:

1. download the latest version of this code.
```
git clone https://github.com/GDMarco/nupropa/
```

2. In the downloaded folder, create a "build/" folder and navigate inside.

3. Install the code with CMake with:
```
cmake ..
make
```

(To note that it requires to adopt both C++14 and C++17 libraries.)

4. If it compiles without any errors, the code is working!

The generated nupropa.py file could be exported in PYTHONPATH.
Alternatively, this .py script has to be specified in your python script.

The hadronisation module requires a functioning version of [PYTHIA 8](https://pythia.org/). (Possibly, it has to be specified in the installation paths.)
To account also for decays of heavy leptons, hadrons or bosons, this code needs to be used in combination with the [CRPYTHIAxDecays](https://github.com/GDMarco/CRPYTHIAxDecays) plug-in.  

## Disclaimer
This plugin is provided by the authors “as is,” without warranty of any kind. Use it with caution and interpret the results carefully. 
For any questions or comments, please, reach out to the authors. 

## Acknowledgements
