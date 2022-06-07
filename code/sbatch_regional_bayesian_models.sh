#!/bin/bash

## Normative modeling
sbatch --array=1-360 /ncf/hcp/data/analyses/myelin/Adversity_Project/code/batch_container_regional_normative_modeling.sh

