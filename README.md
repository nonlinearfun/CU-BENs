# CU-BENs: A structural modeling finite element library
CU-BEN serial version is a nonlinear static and transient dynamic finite element analysis software for structural applications. It supports truss, frame, and discrete Kirchhoff theory (DKT) shell elements. In addtion, CU-BENs also includes built-in, linear acoustic fluid-structure interaction analysis capability.

## Build instructions (Mac OS X)
Please follow the guidelines below for building CU-BENs on Mac OS X.

### Prerequisites:  
- Mac OS X >= 10.10
- Xcode >= 9.0.0
- Suitesparse 5.6.0

### Steps: 
1. Type the following commands in terminal to install necessary packages (Xcode and Suitesparse):
    - Install Xcode:
        ```bash
        xcode-select --install
        ```
    - To install Suitesparse via command line, you would first need to check whether the package Homebrew is installed by running the following: 
        ```bash
        which brew
        ```
    - If Homebrew isn't already installed, run the following command to install Homebrew:"
        ```bash
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
        ```
    - Type the following commands to install suite-sparse:
        ```bash
        brew tap homebrew/science
        brew install suite-sparse
        ```
2. Set up a folder where CU-BENs will be built. In this case, CU-BENs will be installed on your desktop:
    ```bash
     cd Desktop
    ```
    - Clone CU-BENs from Github repository:
    ```bash
     git clone https://github.com/nonlinearfun/CU-BENs.git
    ```
    Note: make sure you have Git installed. You can install Git by running:
    ```bash
    brew install git
    ```
3. Build CU-BENs:
    ```bash
    cd CU-BENs
    Make all
    ```
        
## Build instructions (Ubuntu)
Please follow the guidelines below for building CU-BENs on Ubuntu.

### Prerequisites:  
- Ubuntu >= 14.04
- Suitesparse
    
### Steps: 
1. Run the following commands to update the apt-get cache:
    ```bash
    sudo apt-get update -y
    ```
2. Install Suitesparse:
    ```bash
    sudo apt-get -y libsuitesparse-dev
    ```
3. Set up a folder where CU-BENs will be built. In this case, CU-BENs will be installed on your desktop:
    ```bash
    cd Desktop
    ```
    - Clone CU-BENs from Github repository:
    ```bash
    git clone https://github.com/nonlinearfun/CU-BENs.git
    ```
    Note: make sure you have Git installed. You can install Git by running:
    ```bash
    sudo apt-get install git -y
    ```
4. Build CU-BENs:
    ```bash
    cd CU-BENs
    Make all
    ```
        
A detailed overview of CU-BENs as well as the theory behind the finite element formulation can be found in the tutorial and theory manual, *italicized Introduction to CU-BEN* and *italicized CUBENs theory manual*. Sample input files are provided under the *italicized Sample_Input_Files* directory to exercise the diifferent built-in functions within CU-BENs. Sample input file *italicized model_def_5d_shell.txt* is used to exercise the restart function within CU-BENs. Please be sure you have ran *italicized model_def_5c_shell.txt* in advance. 
