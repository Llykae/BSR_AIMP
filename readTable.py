from numpy import *
import sys
import os

class _Read:
    def __init__(self, method, *args):
        """
        Abstract Class _Read :
        Read something waith the right method. 
            - method : object that is needed to 
                       read.
            - args : arguments for method
        """
        self._method_to_read = method
        self._args_to_read = args

    @property
    def get_result(self):
        return self._method_to_read(*self._args_to_read)

    def __str__(self,):
        pass
        
        

class ReadFile(_Read):
    def __init__(self,fileName=None, directory=None, method_to_read=None):
        """
        Real Class ReadFile
        Method to simplify the reading of files.
            - text : reading results
        """
        self.filename = fileName
        self.directory = directory
        if method_to_read is None:
            method_to_read = open
        if directory[-1] != "/" and fileName[0] != "/":
            directory += "/"
        args_to_read = [directory + fileName, "r"]
        
        super().__init__(method_to_read, *args_to_read)
        self.text = self.get_result.read()

    @property
    def read(self):
        return self.text

    def __str__(self,txt=sys.stdout):
        if self.directory[-1] != "/":
            self.directory += "/"
        return "Reading file from : {}{}".format(self.directory, self.filename)


class Input(_Read):
    def __init__(self,fileName=None,directory=sys.path[0],**properties):
        """
        Real Class Input
        Method to read the properties of calculation
            - fileName : implicite
            - directory : implicite

        """
        if directory[-1] != "/" and fileName[0] != "/":
            directory += "/"
        self.filename = fileName
        self.directory = directory
        self.properties = {}

        if fileName is None:
            assert isinstance(properties,dict)
            self.set_properties(properties)
        elif os.path.exists(directory + fileName):
            self.read_properties()
        else:
            print("\n### Error occurs")
            print(" -> No file found or no properties initialized.\n")


    def read_properties(self):
        pass

    def set_properties(self, **props):
        self.properties.update(props)

    def get_properties(self):
        return self.properties

    def __str__(self):
        if self.directory[-1] != "/":
            self.directory += "/"
        return "Reading Input File from {}{}".format(self.directory, self.filename)


class ReadMatrix(ReadFile):
    def __init__(self,fileName=None,directory=sys.path[0], method_to_read=None):
        """
        Real Class ReadMatrix
        Method to read a Matrix of calculation from DFTB+
            - fileName : implicite
            - directory : implicite

        """
        if directory[-1] != "/" and fileName[0] != "/":
            directory += "/"
        self.filename = fileName
        self.directory = directory

        super().__init__(fileName, directory, method_to_read)
        
    def read_matrix(self):
        pass

    def set_matrix(self, **props):
        pass

    def organize(self):
        return 

    def __str__(self):
        if self.directory[-1] != "/":
            self.directory += "/"
        return "Reading Matrix File from {}{}".format(self.directory, self.filename)



import numpy as np

class readTable :
    def __init__( self, 
                  file_name:str,
                  comment:str = "",
                  verbose:bool = False,
                  replace:float=None,):
        

        self.file_name  = file_name
        self.comment    = comment
        self._verbose   = verbose
        self.data = None
        
        self.__call__(replace)

    def __collect(self, replace=None):
        fichier = open(self.file_name, "r")
        assert fichier
        
        lst = []
        for lines in fichier:
            if lines[0] != "#" :
                tmp = lines.rstrip('\n\r').split()
                for i in range(len(tmp)):
                    if tmp[i] == '':
                        pass
                    else :
                        try:
                            tmp[i] = float(tmp[i])
                        except:
                            pass
                lst.append(tmp)


        fichier.close()

        row_lengths=[]
        for row in lst:
            row_lengths.append(len(row))
        max_length = max(row_lengths)

        for row in lst:
            while len(row) < max_length:
                row.append( replace )

        self.data = lst
        
    def print(self):
    	if self.verbose == True:
    	    print("debug :: Les donnée sont issues du fichier {}"
                  .format(self.file_name))
    	    print("debug :: Commentaire :: {}"
                  .format(self.comment))
    	else:
            print("Need to allow verbose variable")
            
    def extract(self, column_1):
    	return np.array(self.data[:, column_1])

    def __call__(self, replace=None):
    	self.__collect(replace=replace)
    	if self._verbose == True:
            self.print()

    def write( self, 
               newFile:str = "", 
               comment:str = ""):

        M = self.data.shape[1]
        N = len(self.data[0])
        
        fw = open(newFile,"w")
        fw.write("# " + comment + "\n")
        for j in range(N):
            for i in range(M):
                fw.write("{}".format(self.data[i][j]))
                if i < M-1:
                    fw.write("\t")
            fw.write("\n")
        return 0

