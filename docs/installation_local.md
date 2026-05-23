# Local installation (Linux/macOS)

```bash
python -m venv venv && source venv/bin/activate
pip install --upgrade pip

# 1. UASAL config repos
git clone https://github.com/uasal/config_um.git && pip install -e ./config_um
git clone https://github.com/uasal/config_um_wcc.git && pip install -e ./config_um_wcc

# 2. UASAL archive (for stellar spectra)
git clone https://github.com/uasal/uasal_archive.git
export UASAL_ARCHIVE=$PWD/uasal_archive

# 3. budgie
git clone https://github.com/uasal/budgie.git && pip install -e ./budgie

# 4. stp_etc_imaging
git clone https://github.com/uasal/stp_etc_imaging.git
cd stp_etc_imaging && git checkout develop && pip install -e '.[dev,budget]'

# 5. verify
python -c "import stp_etc_imaging, budgie; print(stp_etc_imaging.__version__, budgie.__version__)"
pytest tests/test_target_list.py
```
