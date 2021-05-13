import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing
import pickle
import time
import json
import functools
import datetime
import sys

from trachoma.trachoma_functions import *

def timer(func):
    """Print the runtime of the decorated function"""
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()    # 1
        print(f"-> Running {func.__name__!r}, starting at {datetime.datetime.now()}")
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time    # 3
        print(f"=> Finished {func.__name__!r} in {run_time:.4f} secs\n\n")
        return value
    return wrapper_timer

def loadParameters(BetFilePath, MDAFilePath, PrevFilePath, InfectFilePath, SaveOutput, OutSimFilePath, InSimFilePath, rho, MDA_Cov):

    '''
    Define all required input parameters.

    Parameters:
    -----------
    BetFilePath: str
        This is the path to the input CSV file with the
        random seed and beta to be used for each simulation.

    MDAFilePath: str
        This is the path to the input CSV file with the
        first and last year of the simulations and with
        the first and last year of MDA.

    PrevFilePath: str
        This is the path where the output CSV file with
        the simulated prevalence will be saved.

    InfectFilePath: str
        This is the path where the output CSV file with
        the simulated infection will be saved.

    SaveOutput: bool
        If True, the last state of the simulations will
        be saved in a pickle file. If False, the last
        state of the simulations will not be saved.

    OutSimFilePath: str
        This is the path where the output pickle file with
        the last state of the simulations will be saved. It
        is only required when SaveOutput = True.

    InSimFilePath: str
        This is the path where the input pickle file with
        the last state of the simulations has been saved.
        If this is provided, the code will skip the burnin
        and resume the previous simulations from this state.
        If this is not provided, the code will start new
        simulations from scratch, including the burnin.

    rho: float
        Systematic Adherence

    MDA_Cov: float
        MDA Coverage

    Returns:
    -----------
    sim_params: dict
        simulation parameters.

    params: dict
        general parameters.

    demog: dict
        demography parameters.
    '''

    # load the values of beta and of the random seed
    simparams = pd.read_csv(BetFilePath)
    beta, seed = simparams.iloc[:, 1], simparams.iloc[:, 0]

    # if we are resuming previous simulations we skip
    # the burnin, otherwise we include the burnin
    if InSimFilePath is None:

        burnin = 40 * 52

    else:

        burnin = 0

    # load the simulation years
    timeparams = pd.read_csv(MDAFilePath)
    start_sim_year, end_sim_year = np.int(timeparams.iloc[0, 0]), np.int(timeparams.iloc[0, 1])
    nyears = end_sim_year - start_sim_year + 1

    # output the simulations every 2 months
    sim_start_date = pd.Timestamp(str(start_sim_year) + '-02-01')
    sim_end_date = pd.Timestamp(str(end_sim_year) + '-12-31')
    sim_dates = pd.date_range(start=sim_start_date, end=sim_end_date, freq='2M')
    sim_times = [burnin + np.int(np.round(t)) for t in np.arange(52 / 6, 52 * nyears + 52 / 6, 52 / 6)]
    sim_times = np.array(sim_times[:len(sim_dates)])

    try:

        # try to load MDA dates from JSON-encoded month list in 'mda_vector' in column 5 if present
        # looks like "[202101,202106,202201]"
        mda_vector = json.loads( timeparams.iloc[0, 4] )

        # [ [2021,1], [2021,6], [2022,1] ]
        mda_date_ints = [ [ np.int( str( t )[0:4] ), np.int( str( t )[4:6]  ) ] for t in mda_vector ]

        # [ 2021-01-01 00:00:00, 2021-06-01 00:00:00, 2022-01-01 00:00:00 ]
        mda_dates = [ pd.Timestamp( str( t[0] ) + '-' + str(t[1]).zfill(2) + '-01' ) for t in mda_date_ints ]

        # number of weeks since 2020-01
        start_date_int = [ 2020, 1 ]
        month_addition = 1 if mda_date_ints[0][1] == 6 else 0
        weeks_from_202001 = int( ( ( mda_date_ints[0][0] - start_date_int[0] ) * 52 ) + ( ( ( mda_date_ints[0][1] - start_date_int[1] ) + month_addition ) * ( 52 / 12 ) ) )

        # TODO FIXME account for  equivalent of ( first_mda - start_sim_year )

        # work out the MDA application times from the specified periods - [ 1, 27, 53 ]
        mda_times = np.array( functools.reduce( lambda acc, k: acc + [ burnin + weeks_from_202001 + ( ( k[0] - mda_date_ints[0][0] ) * 52 ) + ( 1 if k[1] == 1 else 26 ) ], mda_date_ints, [] ) )

    except:

        try:
            # load the MDA years
            first_mda, last_mda = np.int(timeparams.iloc[0, 2]), np.int(timeparams.iloc[0, 3])

            # apply the MDA every 12 months
            mda_start_date = pd.Timestamp(str(first_mda) + '-01-01')
            mda_end_date = pd.Timestamp(str(last_mda) + '-01-01')
            mda_dates = pd.date_range(start=mda_start_date, end=mda_end_date, freq='12MS', normalize=True)
            mda_times = [burnin + np.int(np.round(t)) for t in np.arange(1 + 52 * (first_mda - start_sim_year), 52 * nyears + 52, 52)]
            mda_times = np.array(mda_times[:len(mda_dates)])

        except:
            mda_dates = []
            mda_times = []

    print( 'mda_times:', *mda_times )

    # Decide how long simulation you want and when you want MDA to be carried out
    sim_params = dict(
        timesim=burnin + 52 * nyears,      # years total duration of simulation (*52 so in weeks) including burn-in
        burnin=burnin,                     # years burn-in to be discarded (*52 so in weeks)
        Freq_MDA=1,                        # 1 for annual
        N_MDA=len(mda_times),              # no. of rounds of MDA to be carried out in simulation
        MDA_times=mda_times,               # timepoints when MDA occurs
        MDA_dates=mda_dates,               # dates when MDA occurs
        Out_times=sim_times,               # timepoints in output file
        Out_dates=sim_dates,               # dates in output file
        Beta=beta,                         # beta parameter
        Seed=seed,                         # random seed
        n_sim=len(seed)                    # number of simulations
    )

    # General parameters relating to transmission of infection
    params = dict(

        # Population size
        N=1000,

        # Parameters relating to duration of infection period, including ID period
        av_I_duration=2,
        av_ID_duration=200 / 7,
        inf_red=0.45,
        min_ID=11,

        # Parameters relating to duration of disease period
        av_D_duration=300 / 7,
        dis_red=0.3,
        min_D=1,

        # Parameters relating to lambda function - calculating force of infection
        v_1=1,
        v_2=2.6,
        phi=1.4,
        epsilon=0.5,

        # Parameters relating to MDA
        MDA_Cov=MDA_Cov,   # MDA coverage
        MDA_Eff=0.85,  # Efficacy of treatment
        rho=rho        # Correlation parameter for systematic non-compliance function
    )

    # Demography parameters
    demog = dict(
        tau=1 / (40 * 52),  # death rate in weeks^-1
        max_age=60 * 52,    # maximum age in population
        mean_age=20 * 52    # mean age in population
    )

    return sim_params, params, demog

@timer
def Trachoma_Simulation(BetFilePath, MDAFilePath, PrevFilePath, InfectFilePath, SaveOutput=False, OutSimFilePath=None, InSimFilePath=None, rho=0.3, MDA_Cov=0.8):

    '''
    Longitudinal simulations.

    Parameters:
    -----------
    BetFilePath: str
        This is the path to the input CSV file with the
        random seed and beta to be used for each simulation.

    MDAFilePath: str
        This is the path to the input CSV file with the
        first and last year of the simulations and with
        the first and last year of MDA.

    PrevFilePath: str
        This is the path where the output CSV file with
        the simulated prevalence will be saved.

    InfectFilePath: str
        This is the path where the output CSV file with
        the simulated infection will be saved.

    SaveOutput: bool
        If True, the last state of the simulations will
        be saved in a pickle file. If False, the last
        state of the simulations will not be saved.

    OutSimFilePath: str
        This is the path where the output pickle file with
        the last state of the simulations will be saved. It
        is only required when SaveOutput = True.

    InSimFilePath: str
        This is the path where the input pickle file with
        the last state of the simulations has been saved.
        If this is provided, the code will skip the burnin
        and resume the previous simulations from this state.
        If this is not provided, the code will start new
        simulations from scratch, including the burnin.

    rho: float
        Systematic Adherence. Defaults to 0.3.

    MDA_Cov: float
        MDA Coverage. Defaults to 0.8

    Returns:
    -----------
    None.
    '''

    # make sure that the user has provided all the necessary inputs
    if '.csv' not in BetFilePath or '.csv' not in MDAFilePath or '.csv' not in PrevFilePath or '.csv' not in InfectFilePath:

        message = 'Please provide the directory to the CSV files.'

    elif SaveOutput and (OutSimFilePath is None or '.p' not in OutSimFilePath):

        message = 'Please provide the directory to the pickle file, otherwise set SaveOutput = False.'

    else:

        # load all model parameters
        sim_params, params, demog = loadParameters(BetFilePath, MDAFilePath, PrevFilePath, InfectFilePath, SaveOutput,
        OutSimFilePath, InSimFilePath, rho, MDA_Cov)

        if InSimFilePath is None: # start new simulations

            vals = Set_inits(params=params, demog=demog, sim_params=sim_params)  # set initial conditions
            vals = Seed_infection(params=params, vals=vals)  # seed infection

            if sim_params['N_MDA'] != 0:  # create treatment matrix

                previous_rounds = 0  # previous MDA rounds
                Tx_mat = Tx_matrix(params=params, sim_params=sim_params, previous_rounds=previous_rounds)

            else:

                Tx_mat = []

            def multiple_simulations(j):

                out = sim_Ind_MDA(params=params, Tx_mat=Tx_mat, vals=vals, timesim=sim_params['timesim'],
                demog=demog, bet=sim_params['Beta'][j], MDA_times=sim_params['MDA_times'], seed=sim_params['Seed'][j],
                state=None)

                return out

        else: # continue previous simulations

            if sim_params['N_MDA'] != 0:  # create treatment matrix

                previous_rounds = pickle.load(open(InSimFilePath, 'rb'))[0]['N_MDA']  # previous MDA rounds
                Tx_mat = Tx_matrix(params=params, sim_params=sim_params, previous_rounds=previous_rounds)

            else:

                Tx_mat = []

            def multiple_simulations(j):

                vals = pickle.load(open(InSimFilePath, 'rb'))[j]  # load the previous simulations

                out = sim_Ind_MDA(params=params, Tx_mat=Tx_mat, vals=vals, timesim=sim_params['timesim'],
                demog=demog, bet=sim_params['Beta'][j], MDA_times=sim_params['MDA_times'], seed=None,
                state=vals['State'])

                return out

        # run simulations
        num_cores = multiprocessing.cpu_count()

        start_time = time.time()

        out = Parallel(n_jobs=num_cores)(delayed(multiple_simulations)(j) for j in range(sim_params['n_sim']))

        end_time = time.time()

        # save the simulated prevalence in a CSV file
        prevalence_columns = ['Random Generator', 'bet']
        prevalence_columns.extend([d.strftime(format='%m-%Y') for d in sim_params['Out_dates']])
        ddf = pd.DataFrame(columns=prevalence_columns)
        ddf['Random Generator'] = sim_params['Seed']
        ddf['bet'] = sim_params['Beta']

        for i in range(len(out)):

            ddf.iloc[i, 2:] = [out[i]['True_Prev_Disease_children_1_9'][j - 1] for j in sim_params['Out_times']]

        ddf.to_csv(PrevFilePath, index=None)

        # save the simulated infections in a CSV file
        infection_columns = ['Random Generator', 'bet']
        infection_columns.extend([d.strftime(format='%m-%Y') for d in sim_params['Out_dates']])
        idf = pd.DataFrame(columns=infection_columns)
        idf['Random Generator'] = sim_params['Seed']
        idf['bet'] = sim_params['Beta']

        for i in range(len(out)):

            idf.iloc[i, 2:] = [out[i]['True_Infections_Disease_children_1_9'][j - 1] for j in sim_params['Out_times']]

        idf.to_csv(InfectFilePath, index=None)

        # save all the simulated values in a pickle file
        if SaveOutput:

            pickle.dump(out, open(OutSimFilePath, 'wb'))

        message = 'Running time: ' + format((end_time - start_time), '.0f') + ' seconds.'

    print(message)

    return None
