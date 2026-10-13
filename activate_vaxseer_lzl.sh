# Usage: source activate_vaxseer_lzl.sh
_vaxseer_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source /opt/miniconda3/etc/profile.d/conda.sh
export CONDA_ENVS_PATH="${_vaxseer_root}${CONDA_ENVS_PATH:+:${CONDA_ENVS_PATH}}"
conda activate "${_vaxseer_root}/vaxseer_lzl"
# This host has another Conda bin directory ahead of the activated environment.
export PATH="${CONDA_PREFIX}/bin:${PATH}"
hash -r
unset _vaxseer_root
