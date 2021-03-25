# Sample Application for the Contractor for Time-Delayed Signals

## 1. Prerequisites
### 1.1 **Linux and C++11 (or newer)**

### 1.2. **Codac and Ibex**
Follow the [installation instructions for Ibex and Codac](http://codac.io/install/01-installation.html).

### 1.3 **Vibes**
To visualize the application, Vibes is required. You can also install the latest version from the sources available on the [GitHub development repository](https://github.com/ENSTABretagneRobotics/VIBES). For instance, on Linux systems, a fast installation can be made using the following command lines:

```
    sudo apt-get install qt5-default libqt5svg5-dev cmake git
    git clone https://github.com/ENSTABretagneRobotics/VIBES
    cd VIBES/viewer ; mkdir build ; cd build ; cmake .. ; sudo make install
```
You can click on the icon to launch it, or use a terminal. For instance, on Linux systems:

```
    VIBes-viewer
```

## 2. Build the code in this repository
Clone the repository and use cmake to build the application:
```
    git clone https://github.com/Rapha811/ctc_delay_demonstration
    cd ctc_delay_demonstration
    mkdir build -p
    cd build
    cmake ..
    make
```
Afterwards you can launch the application, for example, directly from the terminal:
```
    ./ctc_delay_demonstration
```

## 3. Modify the application
You can modify several properties of the application in the first few lines of the main.cpp file. For example, you can modify the following things:
* ```double dt = 0.2``` Sampling time of the signal (width of the tube slices).
* ```Interval tdomain(-10,10)``` Time domain over which the signal is simulated. Outside this domain the signal is assumed to be zero.
* ```TFunction signal ("(2/(sqrt(3*1.5)*3.14)^(1/4))*(1-(t^2)/(1.5^2))*exp(-(t^2)/(2*(1.5^2)))")``` Simulated signal (Ricker wavelet in this case, but you can use any signal).
* ```IntervalVector sea({{-7.,7.},{-14.,-1.}})``` Dimensions of the environment (x: -7 to 7,y: -14 to -1).
* ```Vector pa{6.5,-6.}, pb{-6,-4.}``` Position of the receiver (pa) and the sound source/emitter (pb).
* ```double velocity = 0.4``` Speed of signal.
* ```double attenuation_coefficient = 0.1``` Attenuation of the signal.

Afterwards, you have to re-compile and launch the application as described before.
