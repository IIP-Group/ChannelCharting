# Simple Channel Charting MATLAB Simulator

Channel charting simulator with principal component analysis (PCA) and Sammon's mapping (SM) for line-of-sight (LoS) and non-LoS scenarios. 
 
(c) 2016-2021 Christoph Studer, Emre Gonultas, and Said Medjkouh

e-mail: studer@ethz.ch, eg566@cornell.edu, sm2685@cornell.edu

### Important information

If you are using this simulator (or parts of it) for a publication, then you must cite the following paper:

Christoph Studer, Said Medjkouh, Emre Gonultas, Tom Goldstein, and Olav Tirkkonen, "Channel charting: Locating users within the radio environment using channel state information," IEEE Access, Vol. 6, pp. 47862-47598, Aug. 2018

_Erratum: The above paper claims to be using the plane-wave model for the vanilla line-of-sight channels; this is incorrect and a spherical wave model was used to generate our simulation results._

In case you decide to use the Sammon's mapping solver for a publication, then you must cite the original FASTA paper:

Tom Goldstein, Christoph Studer, and Richard G. Baraniuk, "A field guide to forward-backward splitting with a FASTA implementation," Technical Report, arXiv preprint: 1411.3406, Nov. 2014; available at https://arxiv.org/pdf/1411.3406.pdf

### How to start a channel charting simulation

Simply run

```sh
channel_charting
```

which starts a simulation with U=2048 user-equipment (UE) locations and a B=32 basestation (BS) massive MIMO antenna. The simulator plots the following system scenario

![](output/scenario.png?raw=true "")

and computes channel charts for PCA and SM. The channel chart for the vanilla LoS model with Sammon's Mapping (SM) should look as follows:

![](output/CC_SM_LoS_U2048_B32.png?raw=true "")

The simulator then computes trustworthiness (TW) and continuity (CT) performance metrics for a pre-specified set of neighborhood sizes and generates the following result:

![](output/TW_CT_SM_LoS_U2048_B32.png?raw=true "")

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Notes

* For the vanilla LoS channel model, you can freely select the system parameters. For the pregenerated Quadriga LoS (QLoS) and non-LoS channels (QNLoS), you can only use the following parameters: U=2048 and B=32. 
* In order to generate your own Quadriga channel vectors, you need to first install the Quadriga software https://quadriga-channel-model.de/ . We do not provide assistance in doing so, but an example script can be found in the `channels` folder. 
* You can easily add your own (or modify our) channel charting functions below the `switch par.DRmethod` statement. 
* You can easily add your own (or modify our) channel-state information (CSI) feature extraction method. 


### Version history

* Version 0.1: studer@ethz.ch - initial version for GitHub release
