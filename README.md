# PILs potential with PiNNAcLe and PiNet2

## Preparing the environment on Alvis2

```bash
ml TensorFlow/2.6.0-foss-2021a-CUDA-11.3.1 ASE/3.22.0-foss-2021a
virtualenv --system-site-packages ~/envs/pinet2-tf26
pip install git+https://github.com/Teoroo-CMC/PiNN.git@pre-release # the branch should be changed after releasing
```
