from . import sweep1d

"""
This set of convenience functions are used to make measuring with the VNA
faster and easier.
"""
# TODO: spec mode settings
# TODO: edits once VNA driver completed
# TODO: docstrings


#########################
# Setup functions
#########################


def resonator_sweep_setup(v1, power=-30, pm_range=200e6, avg=5,
                          bw=1000, npts=2001, centre=7.25e9):
    """
    Function which sets up the VNA with default settings
    for full resonator sweep.

    Args:
        v1 (instrument): VNA instrument to set params on
        power (int): default -30
        pm_range (float): default -30
        avg (int): default 5
        bw (int): default 1000
        npts (int): default 2001
        centre (float): default 7.25e9
    """
    v1.rf_on()
    # v1.spec_mode('off')
    v1.bandwidth(bw)
    v1.avg(avg)
    v1.power(power)
    v1.start(centre - pm_range)
    v1.stop(centre + pm_range)
    v1.npts(npts)


def power_sweep_setup(v1, avg=5, bw=1000, npts=201):
    """
    Function which sets up the VNA with default settings for
    power sweep on a resonator.

    Args:
        v1 (instrument): VNA instrument to set params on
        avg (int): default 5
        bw (int): default 1000
        npts (int): default 201
    """
    v1.rf_on()
    # v1.spec_mode('off')
    v1.bandwidth(bw)
    v1.avg(avg)
    v1.npts(npts)


def gate_sweep_setup(v1, avg=5, bw=1000, npts=201, power=-40):
    """
    Function which sets up the VNA with default settings for power sweep
    on a resonator.

    Args:
        v1 (instrument): VNA instrument to set params on
        avg (int): default 5
        bw (int): default 1000
        npts (int): default 201
        power (int): default -40
    """
    v1.rf_on()
    # v1.spec_mode('off')
    v1.bandwidth(bw)
    v1.avg(avg)
    v1.npts(npts)
    v1.power(power)

#########################
# General functions
#########################


def sweep2d_vna(v1, startf, stopf, stepf,
                sweep_param2, start2, stop2, step2,
                delay=0.01, key="mag", save=True):
    """
    Function which fakes doing a 2d sweep by setting up the vna to
    do a 'hard' frequency sweep over the given range and then executing
    a 'soft' sweep over the other parameter

    Args:
        v1 (instrument): VNA instrument
        startf (float): tarting frequency
        stopf (float): final frequency
        stepf (float): frequency increment
        sweep_param2 (qcodes parameter): second parameter to sweep
        start2 (float): starting value for second parameter
        stop2 (float): final value for second parameter
        step2 (float): step value for second parameter
        delay (float): min time to wait between step of second parameter
        key (str) (default "mag"): string key used to identify what to
            live plot

    Returns:
        dataset, plot
    """
    v1.start(startf)
    v1.stop(stopf)
    npts = int((stopf - startf) / stepf + 1)
    v1.npts(npts)
    dataset, plot = sweep1d(v1.trace, sweep_param2, start2, stop2, step2,
                            delay=delay, key="mag", save=save)
    return dataset, plot


def gates_to_zero(dec_chans):
    """
    Sets all gate channels to 0
    """
    for chan in dec_chans:
        chan(0)
