Instructions for Boinc clients
-----------------------------------

I) Installing the Boinc client
---------------------------------------
1) Downloading the installer
------------------------------
For Windows, download the Boinc+Virtualbox from either http://boinc.berkeley.edu/download.php or https://boinc.berkeley.edu/download_all.php

For Mac and Linux, Boinc and Virtualbox have to be downloaded separately. The Boinc installer can be downloaded from the same links as above. Virtualbox will have to be downloaded from https://www.virtualbox.org/

Alternatively, click on the Green Join button on http://anthgrid.com/dbnupperbound/ and follow the instructions

2) Installing Boinc+Virtualbox
---------------------------------
a) Check that virtualization is enabled on your machine. If it is not enabled, please enable it through the BIOS (this is an unavoidable step since Virtualbox needs this enabled for running properly)
b) Run the installer to install Boinc+Virtualbox on your machine
c) Open the Boinc Manager GUI

As a sanity check, open the file client_state.xml (for eg. in Windows found in the C:/ProgramData/BOINC folder)
and check whether the vm extensions line is indeed <p_vm_extensions_disabled>0</p_vm_extensions_disabled>. 
If it is instead 1, then either virtualization is not enabled yet on your machine, or it was enabled after the Boinc installation. In which case, uninstall Boinc, delete the residual files (for eg. delete the C:/ProgramData/BOINC folder on Windows), reinstall Boinc, and check the client_state.xml file again.

II) Setting up the Boinc client for the project
-----------------------------------------
Within the Boinc manager,
1) Click on Add Project, or Click on Tools > Add Project
2) In the Choose a Project window, enter http://anthgrid.com/dbnupperbound in the Project URL field at the bottom and proceed
3)a) If you are registering for the project, then use the 'No, new user' button and supply your email address and desired password. On proceeding and clicking finish, you will end up at the project website where you can edit your profile.
3)b) If you are already registered with the project (maybe you registered through the website directly, or you are repeating the installation on a second computer, or had to remove the project from the manager and add it back again (which can often be used to fix obscure problems), or some other reason), use the 'Yes, existing user' button, and supply your existing credentials.

III) Setting up computing preferences
-------------------------------------------
Within the GUI manager, go to View > Advanced View 
then Options > Computing Preferences, to open the Preference Settings window
Within the Computing tab
1) Set Usage Limits and When to Suspend settings to appropriate values so the computations do not affect your day to day work
2) In the Other section, it is recommended to set task storing to store 0 days of work and 0 days of additional work, so that your client downloads and processes only 1 job at a time. 
(With non-zero values, there is a possibility that the client will end up downloading too many tasks and putting them in an execution queue, putting a disproportionate responsibility on your machine, and increasing lead times for task completions.

Similarly, detailed preferences can be set in the Network, Disk/Memory and Schedules tab according to your needs.

IV) Executing the first task
--------------------------------------
If there are pending tasks in the project, which can be checked through http://anthgrid.com/dbnupperbound/server_status.php,
within the Advanced view, in the Projects tab, click on the dbnupperbound project, and click on the Update button on the left
(the project Update button can be frequently used for many things, whenever you do not want to wait for the automatic updates generally scheduled once per day)
(the first task may have already started downloading even before clicking on Update)

The required files take a few minutes to download the first time and the download progress can be monitored in the Transfers tab
Next task onwards, Boinc caches the big files, so task downloads take little time and execution starts immediately.

V) Some peculiarities of task completions
---------------------------------------------
The progress bar of a task can be deceptive, since Boinc is just guessing when it will complete. So the progress bar may start rapidly at first, and keep slowing down as it reaches the 100% mark. It is normal for a task to run for anywhere upto 5 hours, but if it is still running beyond that, there is a problem which should be investigated.

VI) Navigating the Project Website
---------------------------------------------
The project website can give a detailed view of the tasks completed by different users, among other features.

Some good-to-have Sections to be added
----------------------------
1) Running things through an account manager like BAM
2) Running the client on headless systems with no manager GUI available
etc.
