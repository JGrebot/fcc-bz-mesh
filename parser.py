
import re
from argparse import ArgumentParser

DESCRIPTION_SCRIPT = """
Formating a .mesh into the two Vertices & Tetrahedra files.
This files reads a mesh with .mesh format and write the two
files needed for the band solver input. Basically this is a parser
of the .mesh mesh file format.

Output files are named: FILENAME_OUT.tetra and FILENAME_OUT.verti
"""

def main_parser_py():
    """
    Formating a .mesh into the two Vertices & Tetrahedra files.

    Basically this is a parser of the .mesh mesh file format.

    See arguments description by typing: python3 parser.py --help
    """
    ##############################################
    # Parsing command lines arguments
    ##############################################
    parser = ArgumentParser(
                description = DESCRIPTION_SCRIPT
                )
    parser.add_argument("-in", 
                        "--filename_in", 
                        dest="filename_in", 
                        type=str, 
                        required=True,
                        help="""Name of .mesh file to transform in UTOX
                        format.""")
    parser.add_argument("-out", 
                        "--filename_out", 
                        dest="filename_out", 
                        type=str,
                        required=True,
                        help="""Name of output files of UTOX format. Two files,
                        one for Vertices, and one for Tetrahedra are
                        written""")
    args = parser.parse_args()

    filename_in = args.filename_in
    filename_out = args.filename_out

    ##############################################
    # Parameters examples
    ##############################################
    # filename_in = "bz.mesh"
    # filename_out = "bz"

    ##############################################
    # Parsing and writing
    ##############################################
    mesh = parse_file(filename_in)
    write_file(filename_out, mesh)
	
	

# set up regular expressions
# use https://regex101.com/ to visualise these if required  
rx_dict = {
    'Vertices': re.compile(r'^ Vertices'),
    'Tetrahedra': re.compile('^ Tetrahedra'),
    'Element': re.compile(r'^ *(?P<P1>[^ ]*) *(?P<P2>[^ ]*) *(?P<P3>[^ ]*) *(?P<P4>[^ ]*)')
}


# Line parser
def _parse_line(line):
    """
    Do a regex search against all defined regexes and return the key and match
    result of the first matching regex

    Never modify this function. Instead just add the corresponding regular
    expression of the marker you aim in the rx_dict dictionnary.

    """
    for key, rx in rx_dict.items():
        match = rx.search(line)
        if match:
            return key, match
    # if there are no matches
    return None, None



# File parser
def parse_file(filepath):
    """
    Parse mesh created by gmsh

    Parameters
    ----------
    filepath : str
        Filepath for file_gmsh to be parsed
    Returns
    -------
    result :    a dict of vertices/tetrahedra

    """
    # Local Variable
    result = {
	'verti': [],
	'tetra': []
    }
    
    # open the file and read through it line by line
    with open(filepath, 'r') as file_gmsh :
        
        line = file_gmsh.readline()
        while line:
            # at each line check for a match with a regex
            key, match = _parse_line(line)

            if key == 'Vertices':
                nb_vertices = int(file_gmsh.readline())
                for ii in range(0,nb_vertices):
                      line_points = file_gmsh.readline()
                      key_verti, match_verti = _parse_line(line_points)
                      result['verti'].append([match_verti.group('P1'), match_verti.group('P2'), match_verti.group('P3')])
			
            if key == 'Tetrahedra':
                nb_tetra = int(file_gmsh.readline())
                for ii in range(0,nb_tetra):
                      line_tetra = file_gmsh.readline()
                      key_tetra, match_tetra = _parse_line(line_tetra)
                      result['tetra'].append([match_tetra.group('P1'), match_tetra.group('P2'), match_tetra.group('P3'), match_tetra.group('P4')])
            line = file_gmsh.readline()
        # END of while
    return result


def write_file(file_out, dict_in):
    f = open(file_out + ".verti", "w")
    for ii in range(0,len(dict_in['verti'])):
         f.write(dict_in['verti'][ii][0] + " " + dict_in['verti'][ii][1]  + " " + dict_in['verti'][ii][2] + "\n")
    f.close()

    f = open(file_out + ".tetra", "w")
    for ii in range(0,len(dict_in['tetra'])):
         f.write(dict_in['tetra'][ii][0] + " " + dict_in['tetra'][ii][1]  + " " + dict_in['tetra'][ii][2] + " " + dict_in['tetra'][ii][3] + "\n")
    f.close()


if __name__ == "__main__":
	main_parser_py()
