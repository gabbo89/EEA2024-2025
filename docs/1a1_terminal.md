---
layout: default
title: Use of the terminal
nav_order: 2
parent: 1. Introduction
description: A comprehensive guide to understanding epigenetics.
published: true
---

The analyses will be performed on a remote server, which can be accessed remotely from any computer with an internet connection. In order to connect to the server a software with SSH (Secure Shell) capabilities is required, that enable the user to interact with the server as if they were physically present. 

For **windows** user's Mobaxterm is the client reccomended, while for **mac** and **linux** users the terminal is the default option.

The ip address and username are already defined (for any issue contact me) 

## **Windows users**

Mobaxterm can be downloaded freely (home edition) from the official website  
[Mobaxterm Download Page](https://mobaxterm.mobatek.net/download.html) 

![moba_download]({{"/assets/images/image.png" | relative_url }})


Once downloaded, the software can be installed following the instructions on screen or on the powerpoint presentation. 

Files can be graphically accessed using the panel on the left side on Mobaxterm. The user can navigate through the directories and select the files they want to use. The files can be copied and pasted to the local computer.


![moba_panel]({{"/assets/images/image-1.png" | relative_url }})


## **Apple users**
To connect to the server and be able to use interactive and graphical software, you need first to have installed on your computer the XQuartz software, which can be downloaded from the official website [XQuartz Download Page](https://www.xquartz.org/).

Next I reccomend to connect using the default terminal application (type *Terminal* on the magnifying glass icon in the top-right corner of the screen) and connect to the server using the following syntax:

`ssh -X {username}@{srv_name}`

`-X` enable the use of interactive devices that use X11 protocol. 

In addition with very new mac devices you may encounter problems with the XQuartz software, in this case you can add the following command line to your `~/.bashrc` file:

`export _JAVA_OPTIONS='-Dsun.java2d.xrender=false -Dsun.java2d.pmoffscreen=false'`

This should solve the problem of black windows (for example on igv) and other graphical issues. If you encounter any problem, remember to comment or remove the section from the .`~/.bashrc` file. 

For any concern, just tell me!

In order to navigate graphically through the directories and select the files, you can use an FTP client like [Cyberduck](https://cyberduck.io/). There are several alternative software options available, such as for example [FileZilla](https://filezilla-project.org/).

You will need to set the server address, the username and the password.