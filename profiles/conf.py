from types import SimpleNamespace
import os

### Set up coef_info for Coef_Manager

wxuas_dir = os.path.join(os.path.expanduser("~"), ".wxuas")


coef_info = SimpleNamespace(USE_AZURE=None, AZURE_CONNECTION_STRING=None, FILE_PATH=None)
# If you do not use Azure, you MUST use a local coefs folder. See oucass.github.io/Profiles/coefs.html for instructions 
coef_info.USE_AZURE="NO"
# If you are using Azure, retrieve the connection strings from your portal.
coef_info.AZURE_CONNECTION_STRING="***REMOVED***"
# If you are NOT using Azure, put the path to the coefs folder here
coef_info.FILE_PATH=wxuas_dir

