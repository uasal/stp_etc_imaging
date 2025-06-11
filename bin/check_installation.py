import os
import importlib.util

def check_module_installations():
    config_stp_spec = importlib.util.find_spec("config_stp")
    config_um_spec = importlib.util.find_spec("config_um")
    config_stp_wcc_spec = importlib.util.find_spec("config_stp_wcc")
    if config_stp_spec is None:
        print("config_stp is not installed.")
    if config_um_spec is None:
        print("config_um is not installed.")
    if config_stp_wcc_spec is None:
        print("config_stp_wcc is not installed.")

    if (config_stp_spec is not None) and (config_um_spec is not None) and (config_stp_wcc_spec is not None):
        print("All configuration modules detected.")

def check_environment_variables():
    if "UASAL_ARCHIVE" in os.environ:
        print("UASAL_ARCHIVE Environment variable detected.")
    else:
        print("UASAL_ARCHIVE not defined.")

if __name__ == "__main__":
    check_module_installations()
    check_environment_variables()