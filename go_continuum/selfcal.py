#!python3
"""Run predefined self-calibration.
"""
from typing import Optional, Sequence
from pathlib import Path
import argparse
import sys

from astropy.table import QTable
from goco_helpers.clean_tasks import image_sn
from go_continuum.data_handler import SelfcalDataManager
from go_continuum.environment import Environ
import astropy.units as u
import goco_helpers.argparse_actions as actions
import goco_helpers.argparse_parents as parents
import numpy as np

def _selfcal_pipe(args: 'argparse.Namespace') -> None:
    """Run the self-calibration pipeline."""
    # Setup environment and read config
    args.log.info('Starting self-calibration')
    environ = Environ(basedir=args.base, check_env=True, continuum=args.base)
    manager = SelfcalDataManager(args.configfile[0], environ, args.log)

    # Concatenate data
    manager.concat_data()
    #vis = manager.concat_uvdata
    # Set config alias
    #config = manager.config

    # Stats
    table = {}

    # Get continuum visibilities
    manager.get_continuum_vis(resume=args.resume)
    manager.init_weights()
    args.log.info('=' * 80)

    # Initial clean
    args.log.info('Initial clean')
    image_info = manager.clean_continuum(nproc=args.nproc[0],
                                         resume=args.resume,
                                         tclean_nsigma=args.tclean_nsigma)
    peak, rms = image_sn(image_info['fitsimage'])
    args.log.info('Image peak: %s', peak)
    args.log.info('Image rms: %s', rms)
    rms_unit = u.uJy/u.beam
    table['iter'] = ['initial']
    table['image'] = [image_info['fitsimage'].name]
    table['threshold'] = (np.array([image_info['thresholds'][-1].value]) *
                          image_info['thresholds'][-1].unit)
    table['peak'] = np.array([peak.value]) * peak.unit
    table['rms'] = np.array([rms.to(rms_unit).value]) * rms_unit
    table['snr'] = [int((peak/rms).to(1).value)]

    # Iterate over solints
    i = 0
    nsigma = 3
    nsigmas = map(lambda x: float(x.strip()),
                  manager.config['selfcal']['threshold_scale'].split(','))
    solints = map(lambda x: x.strip(),
                  manager.config['selfcal']['solint'].split(','))
    for i, (nsigma, solint) in enumerate(zip(nsigmas, solints)):
        args.log.info('=' * 80)
        args.log.info('Iteration: %i', i + 1)
        args.log.info('Cleaning to %fsigma level', nsigma)

        # Clean
        suffix_ending = f'.selfcal{i}'
        image_info = manager.clean_continuum(nproc=args.nproc[0],
                                             nsigma=nsigma,
                                             suffix_ending=suffix_ending,
                                             savemodel='modelcolumn',
                                             tclean_nsigma=args.tclean_nsigma,
                                             resume=args.resume)
        peak, rms = image_sn(image_info['fitsimage'])
        args.log.info('Image %i peak: %s', i, peak)
        args.log.info('Image %i rms: %s', i, rms)
        table['iter'].append(f'{i}')
        table['image'].append(image_info['fitsimage'].name)
        table['threshold'] = np.append(table['threshold'],
                                       image_info['thresholds'][-1])
        table['peak'] = np.append(table['peak'], peak)
        table['rms'] = np.append(table['rms'], rms)
        table['snr'].append(int((peak/rms).to(1).value))

        # Gaincal and applycal
        caltable = image_info['imagename'].with_suffix('.phase.cal')
        args.log.info('Gain table: %s', caltable)
        args.log.info('Solint: %s', solint)
        manager.self_calibrate(caltable, solint, resume=args.resume)

    # Amplitude selfcal?
    if ((ap_solint := manager.config['selfcal'].get('ap', fallback=None))
        is not None):
        # Get model column
        suffix_ending = f'.niter{i + 2}'
        image_info = manager.clean_continuum(nproc=args.nproc[0],
                                             nsigma=nsigma,
                                             suffix_ending=suffix_ending,
                                             tclean_nsigma=args.tclean_nsigma,
                                             savemodel='modelcolumn')
        peak, rms = image_sn(image_info['fitsimage'])
        table['iter'].append(f'{i + 1}')
        table['image'].append(image_info['fitsimage'].name)
        table['threshold'] = np.append(table['threshold'],
                                       image_info['thresholds'][-1])
        table['peak'] = np.append(table['peak'], peak)
        table['rms'] = np.append(table['rms'], rms)
        table['snr'].append(int((peak/rms).to(1).value))

        # Amp selfcal table
        caltable = image_info['imagename'].with_suffix('.amp.cal')
        args.log.info('Gain table: %s', caltable)
        args.log.info('Solint: %s', ap_solint)
        manager.self_calibrate(caltable, ap_solint, calmode='ap',
                               resume=args.resume)

    # Final clean
    suffix_ending = '.final'
    image_info = manager.clean_continuum(nproc=args.nproc[0],
                                         nsigma=nsigma,
                                         tclean_nsigma=args.tclean_nsigma,
                                         suffix_ending=suffix_ending)
    peak, rms = image_sn(image_info['fitsimage'])
    table['iter'].append('final')
    table['image'].append(image_info['fitsimage'].name)
    table['threshold'] = np.append(table['threshold'],
                                   image_info['thresholds'][-1])
    table['peak'] = np.append(table['peak'], peak)
    table['rms'] = np.append(table['rms'], rms)
    table['snr'].append(int((peak/rms).to(1).value))

    # Save table
    table = QTable(table)
    args.log.info('\n%s', table)
    table.write(args.configfile[0].with_suffix('.selfcal.ecsv'),
                format='ascii.ecsv', overwrite=True)

def selfcal(args: Optional[Sequence] = None) -> None:
    """Self-cal main program.

    Args:
      args: Optional. Command line args.
    """
    # Pipe and steps
    pipe = [_selfcal_pipe]

    # Argparse configuration
    args_parents = [parents.logger('debug_selfcal.log')]
    parser = argparse.ArgumentParser(
        description='Automated data self-calibration.',
        add_help=True,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=args_parents,
        conflict_handler='resolve',
    )
    parser.add_argument('--resume', action='store_true',
                        help='Resume unfinished steps')
    parser.add_argument('--tclean_nsigma', action='store_true',
                        help='Use tclean built-in nsigma')
    parser.add_argument('-b', '--base',
                        action=actions.NormalizePath,
                        default=Path('./selfcal'),
                        help='Base directory')
    parser.add_argument('-n', '--nproc', type=int, nargs=1, default=[5],
                        help='Number of processes for parallel steps')
    #parser.add_argument('--skip', nargs='+', choices=list(steps.keys()),
    #                    help='Skip these steps')
    parser.add_argument('--pos', metavar=('X', 'Y'), nargs=2, type=int,
                        help='Position of the representative spectrum')
    #parser.add_argument('--uvdata', action=actions.NormalizePath, nargs='*',
    #                    default=None,
    #                    help='Measurement sets')
    parser.add_argument('configfile', action=actions.CheckFile, nargs=1,
                        help='Configuration file name')
    #parser.set_defaults(manager=None)

    # Read and process
    if args is None:
        args = sys.argv[1:]
    args = parser.parse_args(args)
    for step in pipe:
        step(args)
        args.log.info('=' * 80)

if __name__ == '__main__':
    selfcal(sys.argv[1:])
