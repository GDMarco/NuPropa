# NuPropa

Plugin to the CRPropa code to implement the propagation of highly energetic neutrinos through photon and neutrino backgrounds, e.g. the cosmic microwave and neutrino backgrounds. The production of the secondaries of interaction is treated.

(data/ folder contains the interaction rates employed in the interaction modules. Its position has to be sync in the code.)

## Scientific scenarios



## Installation 

To install *NuPropa*, you will need to have CRPropa 3 installed (at https://github.com/CRPropa/CRPropa3/ the latest version).

To install this plugin:

1. download the latest version of this code.
```
git clone https://github.com/GDMarco/NuPropa/
```

2. In the downloaded folder, create a "build/" folder and navigate inside.

3. Install the code with CMake with:
```
cmake ..
make
```

(To note that it requires to adopt both C++14 and C++17 libraries.)

4. If it compiles without any errors, the code is working!
The generated alpinist.py file could be exported in PYTHONPATH.
Alternatively, this .py script has to be specified in your python script.

## Disclaimer
This code has still to be tested, use it with caution. For any questions or comments, please, reach out the authors. 


## Acknowledgements
