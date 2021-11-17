# COVID-19-simulation

Simulation of early COVID-19 using SIR model and variants (SEIR ...). Made by the Laboratory of Sustainable Life Assessment (GYRO) of the Federal Technologycal University - Parana (UTFPR-ct) in the scope of the project GYRO4Life

# Running the simulation

The code runs based on a csv with the same structure of nc85.csv or oa85.csv files which has a time series of confirmed cases and deaths and metadata information about the region being characterized on the line. Both cases and deaths have to be given for the simulation.

The main code is simulação.py, which receives a couple of arguments:

- 1: region code (for the csv being used). In case the argument is empty ("-"), it will run for all lines of the csv [ex: -28]
- 2: Name of the csv file with confirmed cases (omit the '.csv') [ex: nc85.csv -> -nc85]
- 2: Name of the csv file with confirmed deaths (omit the '.csv') [ex: oa85.csv -> -oa85]
- 3: Fitting method [-0: basinhopp, -1: differential evolution [default], -2: powell, -3: cobyla] [ex: -1]
- 4: Boolean and quantity of opening and closure regimes for the simulation for confirmed cases (works as a contingency method reducing the probability of infection). '-0-0' ignores this factor for a simulation without contingency methods. If a quantity is given on the second argument, the boolean argument must be 1 [ex: '-1-1']
- 5: Boolean and quantity of opening and closure regimes for the simulation for confirmed deaths (works as a contingency method reducing the probability of infection). '-0-0' ignores this factor for a simulation without contingency methods. If a quantity is given on the second argument, the boolean argument must be 1 [ex: '-1-1']
- 6: Type of simulation [-n: simulation of one location (one csv line), -s: simulation of all csv locations, -b: bootstrap of one location [has uncertainty], -sl: simulation of a location with sensibility analysis] [ex: -n]
- 7: Simulation period in days [ex: -200]
- 8: number of days for validation [ex: -5]
- 9: Subtype of simulation [-mod: hospitalization simulation, -std: SEIR simulation with asymptomatic and deaths]
- 10: Run tests and additional graphics [-0: no, -1: yes]

Example call for a SEIR simulation with bootstrap using cases and deaths in Brazil. The simulation is done for 200 days and with a validation of 5 days.

```bash
python simulacao.py -28 -nc85 -oa85 -1 -1-2-0-0 -b -200 -5 -str -0
```
