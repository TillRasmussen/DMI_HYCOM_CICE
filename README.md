# DMI_HYCOM_CICE
This code is the cap for HYCOM, CICE and the driver for the DMI version of HYCOM-CICE
1) Adjust user settings in Make.csh
2) ./Make.csh (automatically calls 'make clean' first.)

Original codes (HYCOM+CICE) are very close to the original ones, however the small changes needs to be checked in order to be included. 
Outstanding issues:
1/ Clean up
2/ There is an issue with domain splittings when HYCOM and CICE delayout do not cover the same domain. I am not sure if this can be handled with tiles

