# node types might be bc, junction, compressor, turbine, etc.
class Coupling_node():
    """ Base class for node dynamics

    Should not e instantiated, serves as a prototype for
    implementaiton of node dynamics.

    Attributes
    ----------
    
    id : int/str
     Unique id of the coupling node
    pde : class
     pde system class, see pde.py for the specifications of
     the required attributes and methods of the pde classes and isoeuler.py
     for an example FV implementation
    nr_pipes : int
     nr of pipes connected to this node
    nus : dict of pipe_id -> np.array each with length 2
     direction of each pipe relative to the node with its 2-norm 
     encodes the thickness of the adjacent pipe.
    dirs : dict of pipe_id -> int
     indication which end of the pipe is connected to the node, if pipe_id:0
     the pipe with id pipe_id is connected with its left-hand side boundary
     to the node and setting pipe_id:1 indicates that the pipe is connected
     with the right-hand side boundary.
    Methods
    -------
    
    empty(id)
     returns an empty boundary object
    
    """
    # cstr
    def __init__(self, id, pde, nus, dirs):
        """
        Parameters
        ----------

        id : int/str
         Unique id of the coupling node
        pde : class
         pde system class, see pde.py for the specifications of
         the required attributes and methods of the pde classes and isoeuler.py
         for an example FV implementation
        nus : dict of pipe_id -> np.array each with length 2
         direction of each pipe relative to the node with its 2-norm 
         encodes the thickness of the adjacent pipe.
        dirs : dict of pipe_id -> int
         indication which end of the pipe is connected to the node, if pipe_id:0
         the pipe with id pipe_id is connected with its left-hand side boundary
         to the node and setting pipe_id:1 indicates that the pipe is connected
         with the right-hand side boundary.
        """
        self.id = id
        self.pde = pde
        self.nr_pipes = len(nus)
        self.nus = nus
        # dirs[pipe.id] = 0 means pipe i is connected with the cells x=0 to the edge
        # and with cells at x=L otherwise
        self.dirs = dirs
    # methods
    @classmethod
    def empty(cls, id):
        cpl_node = cls({}, {})
        cpl_node.id = id
        return cpl_node
    # todo : add calc(pipes, t) base interface
