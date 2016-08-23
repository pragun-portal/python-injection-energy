# python-injection-energy
Python script to calculate the energy in pierce phase and follow through phases of an injection using force-plate data from a csv file.

Usage:
Place file in a folder containing all the csv files that need to be processed. The script will compute the total energy, the energy in the pierce phase and the energy in the follow-through phase for each csv file, and create a csv-file "energy_report.txt" in the current directory. A plot of the jet velocity with the events start of injection, end of injection, start of pierce, end of pierce, start of follow through, and end of follow-through marked. The script excpects the user to do a visual check using the plot and to manually close the plot, after which the script will move to process the next csv file in the folder. The unit used for energy is Joules.

The repository contains two types of python scripts:
- calculate_energy_bento.py
Use this for csv files generated using a Bentobox system.
It uses a moving window of size 100 on account for the different sampling rate on a Bentobox system, and a window size of 10 to calculate the initial offset.
- calculate_energy_oscope.py
Use this for csv files generated using an Oscilloscope
It uses a moving window of size 10, and a window size of 100 to calculate the initial offset.

The following example shows 

```
cd data
python calculate_energy.py
```

The output format of energy_report.txt is as follows:

|CSV Filename | Total energy in injection | Energy in Pierce Phase | Energy in Follow Through Phase |
|-------------|--------------------------:|-----------------------:|-------------------------------:|
|09-01-2015--13-27-44 - WO-100624.csv |	8.12|	2.59|	5.52|

