# Stimulation-induced-decoupling

## Decoupling of interacting neuronal populations by time-shifted stimulation through spike-timing-dependent plasticity

**Abstract**: ‎The synaptic organization of the brain is constantly modified by activity-dependent synaptic plasticity‎. ‎In several neurological disorders‎, ‎abnormal neuronal activity and pathological synaptic connectivity may significantly impair normal brain function‎. ‎Reorganization of neuronal circuits by therapeutic stimulation has the potential to restore normal brain dynamics‎. ‎Increasing evidence suggests that the temporal stimulation pattern crucially determines the long-lasting therapeutic effects of stimulation‎. ‎Here‎, ‎we tested whether a specific pattern of brain stimulation can enable the suppression of pathologically strong inter-population synaptic connectivity through spike-timing-dependent plasticity (STDP)‎. ‎More specifically‎, ‎we tested how introducing a time shift between stimuli delivered to two interacting populations of neurons can effectively decouple them‎. ‎To that end‎, ‎we first used a tractable model‎, ‎i.e.‎, ‎two bidirectionally coupled leaky integrate-and-fire (LIF) neurons‎, ‎to theoretically analyze the optimal range of stimulation frequency and time shift for decoupling‎.  ‎We then extended our results to two reciprocally connected neuronal populations (modules) where inter-population delayed connections were modified by STDP‎. ‎As predicted by the theoretical results‎, ‎appropriately time-shifted stimulation causes a decoupling of the two-module system through STDP‎, ‎i.e.‎, ‎by unlearning pathologically strong synaptic interactions between the two populations‎. ‎Based on the overall topology of the connections‎, ‎the decoupling of the two modules‎, ‎in turn‎, ‎causes a desynchronization of the populations that outlasts the cessation of stimulation‎. ‎Decoupling effects of the time-shifted stimulation can be realized by time-shifted burst stimulation as well as time-shifted continuous simulation‎. ‎Our results provide insight into the further optimization of a variety of multichannel stimulation protocols aiming at a therapeutic reshaping of diseased brain networks‎.

## Usage

- ```two-neuron_motif_LIF.cpp``` reproduces the key results shown in Figure 3 and Figure S3 in the manuscript.
- ```two-module_network_LIF.cpp``` reproduces the key results shown in Figure 4, Figure S1 and Figure S4 in the manuscript.

#### Compilation

- Compile and execute either of the scripts with the following syntax:

```
g++ -std=c++11 two-neuron_motif_LIF.cpp
./a.out

g++ -std=c++11 two-module_network_LIF.cpp -lfftw3
./a.out
```
#### Requirements

- C++11,
- g++ compiler,
- FFTW subroutine library,
- ```random_number.cpp``` and ```spike_statistics.cpp```.

## Citation

- Madadi Asl, M., Valizadeh, A., & Tass, P. A. (2023). Decoupling of interacting neuronal populations by time-shifted stimulation through spike-timing-dependent plasticity. PLOS Computational Biology, 19(2), e1010853. [https://doi.org/10.1371/journal.pcbi.1010853](https://doi.org/10.1371/journal.pcbi.1010853).


