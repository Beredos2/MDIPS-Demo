Massively Distributed In-Parallel Sensing (MDIPS) is the concept of spreading a huge ammount of sensors over a very large volume of space to collect data that can be studied in time independently of space, and viceversa. 
In this example, a small swarm of femotsatellites is set to orbit the Earth in Low Earth Orbit (LEO), and collect data on the magnetic field vector (MFV). After 5 orbits, the simulation stops (the data collected gets huge). 
What follows are some preliminary attempts to try and understand how the data collected can be read, interpreted, and benefited from. 

Programming paradigm: Object oriented

Object classes: celeBod, drone, stereoTaxicSpace, swarm, tiktok. 
-celeBod: a celestial body with some useful properties for calculating some useful astrodynamic quantities, and providing a source for the magnetic field. 
-drone: a femtosatellite class. It is a spacecraft that orbits the Earth and can request the value of the MFV as a function of its position relative to Earth and the aboslute time. 
-stereoTaxicSpace: is an object that represents a virtual grid overlaid on space and that is fixed to the celestial body's frame of reference. 
-swarm: a collection of drones working together to generate the data set. They orbit in a dense group. 
-tiktok: a clock that holds some important time functions and quantities. 

Functions: genData, monitorCel
-genData: ignore, its an incomplete experiment. 
-monitorCel_F: a function that takes as arguments a celestial body, a swarm, a stereotaxic space, and an index that addresses the space, and monitors the addressed cell for all the data set to find and store measurements
               made within it. 
               
