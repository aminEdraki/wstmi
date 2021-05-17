# wSTMI

This is a MATLAB implementation of wSTMI, the speech intelligibility prediction algorithm proposed [here](https://ieeexplore.ieee.org/document/9269417). The algorithm takes the time-aligned clean and degraded/processed speech signals as inputs. The output ```d``` is expected to have a monotonically increasing relation with speech intelligibility.   Please refer to Sec. V of the reference cited below for guidance on interpreting algorithm output.

## Installation
Before using the wSTMI function, the ```utils``` folder has to be added to the search path for the current MATLAB session.


## Usage
The function ```wstmi``` takes three inputs:
```MATLAB
d = wstmi(clean_speech, degraded_speech, sampling_frequency);
```

* ```clean_speech```: An array containing a single-channel clean (reference) speech signal.
* ```degraded_speech```: An array containing a single-channel degraded/processed speech signal.
* ```sampling_frequency```: The sampling frequency of the input signals in ```Hz```.

Note that the clean and degraded speech signals must be time-aligned and of the same length.


## Citing wSTMI
If you use wSTMI, please cite the reference below:
```
@article{edraki2020speech,
  title={Speech Intelligibility Prediction Using Spectro-Temporal Modulation Analysis},
  author={Edraki, Amin and Chan, Wai-Yip and Jensen, Jesper and Fogerty, Daniel},
  journal={IEEE/ACM Transactions on Audio, Speech, and Language Processing},
  volume={29},
  pages={210--225},
  year={2020},
  publisher={IEEE}
}
```

## License
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
