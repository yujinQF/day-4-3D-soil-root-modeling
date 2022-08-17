# RSWMS

v9 is stable.
v10 is in the making.

 ## RSWMS installation
 ### Linux
 ```
 sudo apt-get install update
 sudo apt-get update
 sudo apt-get install build-essential
 sudo apt-get install cmake
 sudo apt-get install gcc
 sudo apt-get install gfortran
 sudo apt-get install openmpi*
 sudo apt-get install liblapack-dev
 sudo apt-get install *mumps*
 ```
 
```
cd ./sparskit
make
cd ..
make
```

### Windows
From the Microsoft Store app, install “Ubuntu 18.04.5 LTS”.

If you have a WSLRegisterDistribution failure then -->
On the Taskbar, click on the Windows Search bar and type Control Panel. Then in the results shown, click on Control Panel.

Now click on Uninstall a Program.
Then, in the left panel of the window, click on Turn Windows Feature on or off

`
! You need admin right to enter this panel
`

Now, scroll down until the end and enable Windows Subsystem for Linux Option.
If the option is already enabled, then disable it and restart your system
 	Then restart your system and then check if the system is clear of the error.
  
In the Ubuntu terminal after installation:

The first time on Ubuntu, you be ask to create a user name (no capital letters and if you use symbols, please use "_") then set a password.
After:

follow the same instruction as for Linux installation.

## Run RSWMS
Follow manual instruction.

