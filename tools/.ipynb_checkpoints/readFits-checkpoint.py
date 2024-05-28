import numpy as np
import os
import astropy.io.fits as fits
from argparse import ArgumentParser

def process_rf_fits(fname, tag, star_cut=.4, rf_bands=None, data_dir='$PSCRATCH/data/run1'):
    if rf_bands is None:
        rf_bands = ['H', 'J', 'Y', 'zpp', 'ip', 'r', 'V', 'B', 'u', 'FLUX_RADIUS', 'PHOTOZ', 'zpec']
    with fits.open(f'{data_dir}/{fname}.fits') as hdul:
        star_class = hdul[1].data['class_star']
        star_filt = star_class < star_cut
        phot_data = np.zeros((np.sum(star_filt), len(rf_bands)+1))
        # phot_err = np.zeros_like(phot_data)
        for i,rb in enumerate(rf_bands):
            phot_data[:,i] = hdul[1].data[rb][star_filt]
            # phot_err[:,i] = hdul[1].data[rb+'err'][star_filt]
    phot_data[:,-1] = tag
    # phot_err[:,-1] = tag
    header = rf_bands+['blend']
    return phot_data, header

def process_rf_vali(fname, bflag, star_cut=.4, rf_bands=None, data_dir='$PSCRATCH/data/run1'):
    if rf_bands is None:
        rf_bands = ['H', 'J', 'Y', 'zpp', 'ip', 'r', 'V', 'B', 'u', 'FLUX_RADIUS', bflag, 'bmu_dist', 'PHOTOZ', 'zpec']
    with fits.open(f'{data_dir}/{fname}.fits') as hdul:
        star_class = hdul[1].data['class_star']
        star_filt = star_class < star_cut
        phot_data = np.zeros((np.sum(star_filt), len(rf_bands)))
        phot_err = np.zeros_like(phot_data)
        for i,rb in enumerate(rf_bands):
            if rb == bflag:
                phot_data[:,i] = (hdul[1].data[rb][star_filt] > 0).astype(int)
            else:
                phot_data[:,i] = hdul[1].data[rb][star_filt]
            # phot_err[:,i] = hdul[1].data[rb+'err'][star_filt]
    header = rf_bands[:-1]+['blend']
    return phot_data, header

def process_rf_vali_err(fname, bflag, star_cut=.4, rf_bands=None, data_dir='$PSCRATCH/data/run1'):
    if rf_bands is None:
        bands = ['H', 'J', 'Y', 'zpp', 'ip', 'r', 'V', 'B', 'u']
        berr = [b+'ERR' for b in bands]
        rf_bands = bands + berr +  ['FLUX_RADIUS', bflag, 'bmu_dist']
    with fits.open(f'{data_dir}/{fname}.fits') as hdul:
        star_class = hdul[1].data['class_star']
        star_filt = star_class < star_cut
        phot_data = np.zeros((np.sum(star_filt), len(rf_bands)))
        phot_err = np.zeros_like(phot_data)
        for i,rb in enumerate(rf_bands):
            if rb == bflag:
                phot_data[:,i] = (hdul[1].data[rb][star_filt] > 0).astype(int)
            else:
                phot_data[:,i] = hdul[1].data[rb][star_filt]
            # phot_err[:,i] = hdul[1].data[rb+'err'][star_filt]
    header = rf_bands[:-1] + ['blend']
    return phot_data, header

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-n", type=int)
    args = parser.parse_args()
    pscratch = os.environ['PSCRATCH']
    data_dir = f'{pscratch}/data/run{args.n}'
    out_dir = f'{pscratch}/output/run{args.n}'

#     vali_strong_err, err_header = process_rf_vali_err('vali_strong', 'b1', star_cut=1, data_dir=data_dir)
#     csv_headstr = ','.join(err_header)
#     np.savetxt(f'{out_dir}/vali_strong_err_messy.csv', vali_strong_err, delimiter=',', header=csv_headstr, comments='')

    # training_pure, _ = process_rf_fits('train_pure', tag=0, star_cut=1, data_dir=data_dir)
    # training_blend, _ = process_rf_fits('iden_weak_unrec', tag=1, star_cut=1, data_dir=data_dir)  
    # training_strong, _ = process_rf_fits('iden_strong_unrec', tag=1, star_cut=1, data_dir = data_dir)
    # vali_weak, _ = process_rf_vali('vali_mapped', 'b2', star_cut=1, data_dir=data_dir)
    # vali_strong, _ = process_rf_vali('vali_mapped', 'b1', star_cut=1, data_dir=data_dir)

    # csv_header_vali =  ['H', 'J', 'Y', 'zpp', 'ip', 'r', 'V', 'B', 'u', 'FLUX_RADIUS', 'blend', 'bmu_dist', 'pz', 'sz']
    # csv_headstr_vali = ','.join(csv_header_vali)
    # 
    # csv_header =  ['H', 'J', 'Y', 'zpp', 'ip', 'r', 'V', 'B', 'u', 'FLUX_RADIUS', 'pz', 'sz', 'blend']
    # csv_headstr = ','.join(csv_header)
    # 
    # np.savetxt(f'{out_dir}/training_pure_messy.csv', training_pure, delimiter=',', header=csv_headstr, comments='')
    # np.savetxt(f'{out_dir}/training_blend_messy.csv', training_blend, delimiter=',', header=csv_headstr, comments='')
    # np.savetxt(f'{out_dir}/training_strong_messy.csv', training_strong, delimiter=',', header=csv_headstr, comments='')
    # np.savetxt(f'{out_dir}/vali_weak_som.csv', vali_weak, delimiter=',', header=csv_headstr_vali, comments='')
    # np.savetxt(f'{out_dir}/vali_strong_som.csv', vali_strong, delimiter=',', header=csv_headstr_vali, comments='')

    rf_bands_weak = ['H', 'J', 'Y', 'zpp', 'ip', 'r', 'V', 'B', 'u', 'FLUX_RADIUS', 'b2', 'bmu_dist', 'PHOTOZ', 'zpec', 'somz']
    rf_bands_strong = ['H', 'J', 'Y', 'zpp', 'ip', 'r', 'V', 'B', 'u', 'FLUX_RADIUS', 'b1', 'bmu_dist', 'PHOTOZ', 'zpec', 'somz']

    vali_weak, _ = process_rf_vali('vali_mapped_somz', 'b2', star_cut=1, rf_bands=rf_bands_weak, data_dir=data_dir)
    vali_strong, _ = process_rf_vali('vali_mapped_somz', 'b1', star_cut=1, rf_bands=rf_bands_weak, data_dir=data_dir)

    csv_header_vali =  ['H', 'J', 'Y', 'zpp', 'ip', 'r', 'V', 'B', 'u', 'FLUX_RADIUS', 'blend', 'bmu_dist', 'pz', 'sz', 'somz']
    csv_headstr_vali = ','.join(csv_header_vali)
    
    np.savetxt(f'{out_dir}/vali_weak_somz.csv', vali_weak, delimiter=',', header=csv_headstr_vali, comments='')
    np.savetxt(f'{out_dir}/vali_strong_somz.csv', vali_strong, delimiter=',', header=csv_headstr_vali, comments='')

    
