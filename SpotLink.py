import sys
import multiprocessing
from System.MSData import CDataPack
from Function.MSFunctionIO import CFunctionConfigParser
from Task.MSTaskIO import CTaskPreparing, CTaskPrepareIni, CTaskPrePareMS2


class CStaff:

    def __init__(self, input_args):
        self.data_pack = CDataPack()
        self.argv = input_args

    def start(self):
        self.__captain_run_flow()
        self.__prepareUserInputs(data_package)

    def __captain_run_flow(self):
        function_config = CFunctionConfigParser()
        function_config.file_to_config(self.argv[1], self.data_pack.my_config)

    def __prepareUserInputs(data_package):
        prepare_task = CTaskPreparing()
        prepare_task.work(data_package)
        ini_task = CTaskPrepareIni()
        ini_task.work(data_package)
        ms2_task = CTaskPrePareMS2()
        ms2_task.work(data_package)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    staff = CStaff(sys.argv)
    staff.start()