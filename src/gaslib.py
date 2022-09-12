from network import *
from zipfile import ZipFile
import xmlschema
import tempfile
import pint


def gaslib_reader(input_network_zip, input_schema_zip):
    tmp_target_path = tempfile.mkdtemp()
    # read xsd schema
    with ZipFile(input_schema_zip,"r") as zip_ref:
        zip_ref.extractall(tmp_target_path)
    schema = xmlschema.XMLSchema(tmp_target_path+'/schemaDefinition/'+'Gas.xsd')
    # read network topology
    with ZipFile(input_network_zip,"r") as zip_ref:
        zip_ref.extractall(tmp_target_path+'/nw')
    # todo: add exception if false: schema.is_valid(tmp_target_path+'/nw/GasLib-11.net'))
    return schema.to_dict(tmp_target_path+'/nw/GasLib-11.net')


# def construct_edges_isoeuler(nodes_dict, edges_dict):
#     nodes = {}
#     edges = {}
#     res_edges += construct_pipe_edges(edges_dict['pipe'], edges_dict)
#     res = construct_valve_edges(edges_dict['valve'], edges_dict)
#     res += construct_compressore_nodes(edges_dict['compressorStation'], edges_dict)
#     # add more node types if necessary
#     return res
#     for key in nodes_dict:
#         if (key == 'pipe'):
#             pass
#         elif (key == 'valve'):
#             pass
#         elif (key == 'compressorStation'):
#             pass
#         # add more edge types here if necessary
#         else:
#             #error
#             pass
            

# def cstr_edges_isoeuler(edge_type):
#     nodes = {}
#     edges = {}  
#     if (edge_type == 'pipe'):
#         return cstr_pipe_isoeuler
#     elif (edge_type == 'valve'):
#     elif (edge_type == 'compressorStation'):
#     else:
#         # todo exception
#         pass

# def update_nodes(nodes, edges, all_nodes):
#     for n in nodes:
#         if all

# todo everything with a suffix "isoeuler" goes into a gaslib reader class

# todo postprocessing setting functions for source, sinks, compressors, valves
# todo entries in the middle of a pipe are replace by junctions
# todo: move physics of the isoeuler in separate class (i.e separate from cweno)
ureg_isoeuler = pint.UnitRegistry()

dx_isoeuler = 100 * ureg_isoeuler.meter
cfl_isoeuler = 0.7
alpha_isoeuler = 1 
gamma_isoeuler = 1
init_table_isoeuler = {'default' : lambda x : [1.0, 0.0]} # pipe_id -> function, default -< function


# transforms pint unit object to meter floats 
def meter_float(x):
    return float(x.to_base_units().m)

def construct_boundary_isoeuler(gaslib_data_dict, nu, dir):    
    return Boundary(gaslib_data_dict['@id'], [nu], [dir])

def construct_junction_isoeuler(gaslib_data_dict, nu, dir):
    return Junction(gaslib_data_dict['@id'], [nu], [dir])

def construct_nonleaf_sinksource(gaslib_data_dict, nu, dir):
    return Junction(gaslib_data_dict['@id'], [nu], [dir])

def construct_pipe_isoeuler(gaslib_data_dict):
    id = gaslib_data_dict['@id']
    L = gaslib_data_dict['length']['@value'] * ureg_isoeuler(gaslib_data_dict['length']['@unit'])
    if id in init_table_isoeuler:
        init = init_table_isoeuler[id]
    else:
        init = init_table_isoeuler['default']
    return {}, {id : Pipe(id, cweno3(id,
                               math.floor(meter_float(L/dx_isoeuler)),
                               meter_float(L),
                               cfl_isoeuler,
                               alpha_isoeuler,
                               gamma_isoeuler,
                               init))}


def construct_valve_isoeuler(gaslib_data_dict):
    pass

def construct_compressor_isoeuler(gaslib_data_dict):
    pass


gaslib_dict_isoeuler = gaslib_reader('/home/aleks/Downloads/GasLib-11.zip', '/home/aleks/Downloads/schemaDefinition.zip')



node_table_isoeuler = {'sink' : construct_boundary_isoeuler,
                       'source' : construct_boundary_isoeuler,
                       'innode' : construct_junction_isoeuler,
                       'inSinkSource' : construct_nonleaf_sinksource }

edge_table_isoeuler = {'pipe' : construct_pipe_isoeuler,
                       'valve' : construct_valve_isoeuler,
                       'compressorStation' : construct_compressor_isoeuler}

# these pipe diameters are artifical! and have to be set to some value, in this case one meter
valve_diam_isoeuler = 1
compressor_diam_isoeuler = 1


def search_node_isoeuler(node_id):
    for key, nodes_list in gaslib_dict_isoeuler['framework:nodes'].items():
        filtered_nodes = list(filter(lambda n: n['@id'] == node_id, nodes_list))
        if filtered_nodes != []:
            return key, filtered_nodes[0]

def get_edge_length(n1_dict, n2_dict):
    x = np.array([n2_dict['@x'] - n1_dict['@x'], n2_dict['@y'] - n1_dict['@y']], dtype='d')
    return  np.linalg.norm(x)

        
def calc_nus(n1_dict, n2_dict, diam):
    nu1 = np.array([n2_dict['@x'] - n1_dict['@x'], n2_dict['@y'] - n1_dict['@y']], dtype='d')
    nu1 = nu1/get_edge_length(n1_dict, n2_dict)*diam
    nu2 = -nu1
    return [nu1, nu2]

def add_node(all_nodes, node_type, node_dict, nu, dir):
    id = node_dict['@id']
    if id in all_nodes:
        # replace an (existant) source or sink in the middle of a pipe by a junction
        if node_type == 'sink' or node_type == 'source':
            all_nodes[id] = node_table_isoeuler['inSinkSource'](node_dict, nu, dir)
        all_nodes[id].nr_pipes += 1
        all_nodes[id].nus.append(nu)
        all_nodes[id].dirs.append(dir)        
    else:
        all_nodes[id] = node_table_isoeuler[node_type](node_dict, nu, dir)

def add_edge(all_edges, all_nodes, edge_type, edge_dict):
    nodes, edges = edge_table_isoeuler[edge_type](edge_dict)
    # merging the new constructed dicts with dicts containg all
    all_nodes = {**all_nodes, **nodes}
    all_edges = {**all_edges, **edges}

def get_diam_isoeuler(edge_type, edge):
    diam = 0
    if (edge_type == 'valve'):
        diam = valve_diam_isoeuler
    elif (edge_type == 'compressorStation'):
        diam = compressor_diam_isoeuler
    elif (edge_type == 'pipe'):
        diam = edge['diameter']['@value'] * ureg_isoeuler(edge['diameter']['@unit'])
        diam = meter_float(diam)
    return diam



# this method is a general one however the specific child class (e.g. isoeuler)
# has to provide the node_table!
def add_edges_with_nodes(edge_type, edge, all_nodes, all_edges):
    n1_type, n1_dict  = search_node_isoeuler(edge['@from'])
    n2_type, n2_dict  = search_node_isoeuler(edge['@to'])
    diam = get_diam_isoeuler(edge_type, edge)
    L = get_edge_length(n1_dict, n2_dict)
    nus = calc_nus(n1_dict, n2_dict, diam)
    
    add_node(all_nodes, n1_type, n1_dict, nus[0], 0)
    add_node(all_nodes, n2_type, n2_dict, nus[1], 1)

    # an edge might be divided into two, hence passing the all_nodes dict
    add_edge(all_edges, all_nodes, edge_type, edge)

    # set cur nodes and edges
    # update all_nodes, if exists replace (necessary for junctions since we proceed edgewise)
    

    
    
def gaslib_to_isoeuler():
    all_nodes = {}
    all_edges = {}
            
    for edge_type, edges in gaslib_dict_isoeuler['framework:connections'].items():
        if (isinstance(edges, list)):
            for edge in edges:
                add_edges_with_nodes(edge_type, edge, all_nodes, all_edges)
        else: # edges might be a single pipe that is not a list
            add_edges_with_nodes(edge_type, edges, all_nodes, all_edges)

    return all_nodes, all_edges


nodes, edges = gaslib_to_isoeuler()
