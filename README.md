# IBCsignature

### OPTION 1; IBC signature standalone GUI (no need to install MATLAB)
Reproduces the analysis reported in Zare *et al*. To score novel IBC/nonIBC dataset, see option 2 (scoring novel dataset in this GUI will be incorporated in the near future)
#### Supported Operating Systems
###### Windows
Has been successfully tested on Windows 7 and Windows 10 (see installation process below)
###### Mac
See option2 below
###### Linux
See option2 below

#### Installation
-Download *'IBCsignatureGUI.exe'* and *'Data.tar.gz'* file (contains data used in Zare *et al*). If you don't have a software to uncompress '.tar' and '.gz' files, download *'UnZipUnTAr.exe'*.  
-Install *'IBCsignatureGUI.exe'* file as you would any other windows executable file. During the installation process, there will be a one time request to install MATLAB Compiler if not found in your system (MATLAB Runtime Version 9.5 for APPs compiled in MATLAB R2018b). After succesfully installing *'IBCsignatureGUI.exe'*, if need be, install *'UnZipUnTar.exe'*. 
#### Analysis
-Start the installed *'IBCsignatureGUI.exe'* software (the graphical user interface (GUI) is shown below. Note, the cyan *'Initiate G59 analysis & random forest'* button will be invisible and will only appear when data to be analyzed is selected using the knob).  
-To uncompress *'Data.tar.gz'*, start *'UnZipUnTar.exe'* GUI, click *'Uncompress .tar.gz file'* and follow instructions. 
-Click *'Data directory'* button and select uncompressed Data directory/folder.   
-Turn the knob to select data to analyze. Click *'Initiate G59 analysis & random forest'* button to initiate analysis

 
![image](https://user-images.githubusercontent.com/68044059/128934089-49080c28-2775-40e6-b32e-f4e2091f044e.png)


### OPTION 2; Run directly in MATLAB
-Technically, should work on any Operating System with MATLAB installed.
#### Dependency
-Statistics and Machine Learning Toolbox, version 11.4 or higher
-Bioinformatics Toolbox, version 4.11 or higher
-MatSurv (Optional) https://github.com/aebergl/MatSurv
#### Installation and analysis
-Download *
