import os
import contextlib
import pandas as pd
from profiles.conf import coef_info
from abc import abstractmethod


class Coef_Manager_Base:
    @abstractmethod
    def get_tail_n(self, copterID):
        """ Get the tail number corresponding to a short ID number.

        :param int copterID: the short ID number of the copter
        :rtype: str
        :return: the tail number
        """
        pass

    @abstractmethod
    def get_sensors(self, scoopID):
        """ Get the sensor serial numbers for the given scoop.

        :param str scoopID: The scoop's identifier
        :rtype: dict
        :return: sensor numbers as {"imet1":"", "imet2":"", "imet3":"", "imet4":"",\
                                    "rh1":"", "rh2":"", "rh3":"", "rh4":""}
        """
        pass

    @abstractmethod
    def get_coefs(self, type, serial_number):
        """ Get the coefs for the sensor with the given type and serial number.

        :param str type: "Imet" or "RH" or "Wind"
        :param str serial_number: the sensor's serial number
        :rtype: dict
        :return: information about the sensor, including offset OR coefs and calibration equation
        """
        pass


class Coef_Manager(Coef_Manager_Base):
    """ Reads the .profilesrc file to determine if the coefs are in the \
       local file system or on Azure and to determine the file path \
       or connection string, then ingests the data from the proper source. \
       This object can then be queried by scoop number (to get sensor numbers), \
       by sensor numbers (to get coefs), or by copter number (to get tail \
       number).
    """

    def __init__(self):
        """ Create Coef_Manager
        """
        # The sub_manager will implement all abstract methods
        self.sub_manager = None
        if coef_info.USE_AZURE.upper() in "YES":
            pass
            # try:
            #     with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
            #         table_service = TableServiceClient.from_connection_string(conn_str=coef_info.AZURE_CONNECTION_STRING)
            # except Exception as e:
            #     print(e)
            #     raise Exception("There are no valid connection strings in __init__.py")
            # # If this point has been reached, we have an active table_service to pull data from
            # self.sub_manager = Azure_Coef_Manager(table_service)
        elif coef_info.USE_AZURE.upper() in "NO":
            if len(coef_info.FILE_PATH) > 0:
                if os.path.exists(coef_info.FILE_PATH):
                    self.sub_manager = CSV_Coef_Manager(coef_info.FILE_PATH)
                else:
                    raise Exception("The FILE_PATH in conf.py is not valid")
            else:
                raise Exception("When USE_AZURE is NO, the FILE_PATH must be set in conf.py.")
        else:
            raise Exception("USE_AZURE in __init__.py must be YES or NO")

    def get_tail_n(self, copterID):
        """ Get the tail number corresponding to a short ID number.

        :param int copterID: the short ID number of the copter
        :rtype: str
        :return: the tail number
        """
        return self.sub_manager.get_tail_n(copterID)

    def get_sensors(self, scoopID):
        """ Get the sensor serial numbers for the given scoop.

        :param str scoopID: The scoop's identifier
        :rtype: dict
        :return: sensor numbers as {"imet1":"", "imet2":"", "imet3":"", "imet4":"",\
                                    "rh1":"", "rh2":"", "rh3":"", "rh4":""}
        """
        return self.sub_manager.get_sensors(scoopID)

    def get_coefs(self, type, serial_number):
        """ Get the coefs for the sensor with the given type and serial number.

        :param str type: "Imet" or "RH" or "Wind"
        :param [str or int] serial_number: the sensor's serial number
        :rtype: dict
        :return: information about the sensor, including offset OR coefs and calibration equation
        """
        return self.sub_manager.get_coefs(type, str(serial_number))


class CSV_Coef_Manager(Coef_Manager_Base):

    def __init__(self, file_path):
        self.file_path = file_path
        self.coefs = pd.read_csv(os.path.join(file_path, 'MasterCoefList.csv'))
        self.copternums = pd.read_csv(os.path.join(file_path, 'copterID.csv'), names=['id', 'tail'])

    def get_tail_n(self, copterID):
        """ Get the tail number corresponding to a short ID number.

        :param int copterID: the short ID number of the copter
        :rtype: str
        :return: the tail number
        """

        return self.copternums['tail'][self.copternums['id'] == str(int(copterID))].values[0]

    def get_sensors(self, scoopID):
        """ Get the sensor serial numbers for the given scoop.

        :param str scoopID: The scoop's identifier
        :rtype: dict
        :return: sensor numbers as {"imet1":"", "imet2":"", "imet3":"", "imet4":"",\
                                    "rh1":"", "rh2":"", "rh3":"", "rh4":""}
        """
        scoop_info = pd.read_csv(os.path.join(self.file_path, 'scoops.csv'))

        # max_date="0000-00-00"
        # dates = scoop_info.validFrom
        # for date in dates:
        #     if date > max_date:
        #         max_date = date
        #
        # # most_recent = scoop_info[scoop_info.validFrom == max_date]
        # return {"imet1":str(most_recent.imet1.values[0]), "imet2":str(most_recent.imet2.values[0]),
        #         "imet3":str(most_recent.imet3.values[0]), "imet4":None, "rh1":str(most_recent.rh1.values[0]),
        #         "rh2":str(most_recent.rh2.values[0]),  "rh3":str(most_recent.rh3.values[0]), "rh4":None}

        sensors = scoop_info[scoop_info.name == scoopID]

        return {"imet1": str(sensors.imet1.values[0]), "imet2": str(sensors.imet2.values[0]),
                "imet3": str(sensors.imet3.values[0]), "imet4": None, "rh1": str(sensors.rh1.values[0]),
                "rh2": str(sensors.rh2.values[0]), "rh3": str(sensors.rh3.values[0]), "rh4": None}

    def get_coefs(self, stype, serial_number):
        """ Get the coefs for the sensor with the given type and serial number.

        :param str stype: "Imet" or "RH" or "Wind"
        :param [str or int] serial_number: the sensor's serial number
        :rtype: dict
        :return: information about the sensor, including offset OR coefs and calibration equation
        """
        try:
            serial_number = int(serial_number)
        except ValueError:
            # This means it was a tail number, which can't be cast  to int
            pass

        coefs = self.coefs
        a = coefs.A[coefs.SerialNumber == str(serial_number)][coefs.SensorType == stype].values[0]
        b = coefs.B[coefs.SerialNumber == str(serial_number)][coefs.SensorType == stype].values[0]
        c = coefs.C[coefs.SerialNumber == str(serial_number)][coefs.SensorType == stype].values[0]
        d = coefs.D[coefs.SerialNumber == str(serial_number)][coefs.SensorType == stype].values[0]
        eq = coefs.Equation[coefs.SerialNumber == str(serial_number)][coefs.SensorType == stype].values[0]
        offset = coefs.Offset[coefs.SerialNumber == str(serial_number)][coefs.SensorType == stype].values[0]

        return {"A":a, "B":b, "C":c, "D":d, "Equation":eq, "Offset":offset}
