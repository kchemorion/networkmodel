import numpy as np
import pandas as pd
import libsbml

def sanitize_id(name):
    """ Sanitize identifiers to be valid SBML IDs """
    return name.replace(' ', '_').replace('-', '_').replace('/', '_').replace('β', 'beta').replace('γ', 'gamma').replace('α', 'alpha')

def create_matrices(filename):
    """ Create activation and inhibition matrices from an Excel file """
    nw = pd.read_excel(filename)
    num_of_nodes = len(nw['Nodes'])
    mact = np.zeros((num_of_nodes, num_of_nodes))
    minh = np.zeros((num_of_nodes, num_of_nodes))
    node_names = nw['Nodes']
    for i in range(num_of_nodes):
        act_in_line = [sanitize_id(node.strip()) for node in str(nw['Activators'][i]).split(',')]
        inh_in_line = [sanitize_id(node.strip()) for node in str(nw['Inhibitors'][i]).split(',')]
        mact[i, :] = [1 if node in act_in_line else 0 for node in node_names]
        minh[i, :] = [1 if node in inh_in_line else 0 for node in node_names]
    return mact, minh, node_names.tolist(), num_of_nodes

def add_plot_annotations(model, groups, node_names):
    """ Add plot annotations to the SBML model """
    plot_annotation = libsbml.XMLNode(libsbml.XMLTriple("plot_config", "", ""), libsbml.XMLAttributes())
    for group_name, indices in groups.items():
        group_node = libsbml.XMLNode(libsbml.XMLTriple("plot_group", "", ""), libsbml.XMLAttributes())
        group_node.addAttr("name", group_name)
        group_node.addAttr("type", "timeSeries")
        for idx in indices:
            species_id = sanitize_id(node_names[idx])
            species_node = libsbml.XMLNode(libsbml.XMLTriple("species", "", ""), libsbml.XMLAttributes())
            species_node.addAttr("id", species_id)
            group_node.addChild(species_node)
        plot_annotation.addChild(group_node)
    if model.isSetAnnotation():
        existing_annotation = model.getAnnotation()
        existing_annotation.addChild(plot_annotation)
    else:
        model_annotation = libsbml.XMLNode()
        model_annotation.addChild(plot_annotation)
        model.setAnnotation(model_annotation)

def export_to_sbml_with_params(filename, node_names, mact, minh, xinit, clamped, gamma, h):
    sbml_document = libsbml.SBMLDocument(3, 1)
    model = sbml_document.createModel()
    compartment = model.createCompartment()
    compartment.setId('default')
    compartment.setConstant(True)
    compartment.setSize(1.0)

    decay_param = model.createParameter()
    decay_param.setId('gamma')
    decay_param.setValue(gamma)
    decay_param.setConstant(True)

    steepness_param = model.createParameter()
    steepness_param.setId('h')
    steepness_param.setValue(h)
    steepness_param.setConstant(True)

    sanitized_node_ids = [sanitize_id(node) for node in node_names]
    for i, node in enumerate(sanitized_node_ids):
        species = model.createSpecies()
        species.setId(node)
        species.setCompartment('default')
        species.setInitialConcentration(xinit[i])
        species.setBoundaryCondition(bool(clamped[i]))
        species.setConstant(False)
        species.setHasOnlySubstanceUnits(False) 

    # Add reactions
    for i, target_node in enumerate(sanitized_node_ids):
        for j, activator_node in enumerate(sanitized_node_ids):
            if mact[j, i]:
                reaction = model.createReaction()
                reaction.setId(f"{activator_node}_activates_{target_node}")
                reaction.setReversible(False)
                reaction.setFast(False)
                reactant = reaction.createReactant()
                reactant.setSpecies(activator_node)
                reactant.setStoichiometry(1)
                reactant.setConstant(False)  # Adding the 'constant' attribute
                product = reaction.createProduct()
                product.setSpecies(target_node)
                product.setStoichiometry(1)
                product.setConstant(False)  # Adding the 'constant' attribute
                kinetic_law = reaction.createKineticLaw()
                formula = f"((-exp(0.5 * h) + exp(-h * ({activator_node} - 0.5))) / ((1 - exp(0.5 * h)) * (1 + exp(-h * ({activator_node} - 0.5))))) - gamma * {target_node}"
                math_ast = libsbml.parseL3Formula(formula)
                kinetic_law.setMath(math_ast)
                                
        for j, inhibitor_node in enumerate(sanitized_node_ids):
            if minh[j, i]:
                reaction = model.createReaction()
                reaction.setId(f"{inhibitor_node}_inhibits_{target_node}")
                reaction.setReversible(False)
                reaction.setFast(False)  # Set the 'fast' attribute to False
                reactant = reaction.createReactant()
                reactant.setSpecies(inhibitor_node)
                reactant.setStoichiometry(1)
                reactant.setConstant(True) 
                product = reaction.createProduct()
                product.setSpecies(target_node)
                product.setStoichiometry(-1)
                product.setConstant(True) 
                kinetic_law = reaction.createKineticLaw()
                formula = f"((-exp(0.5 * h) + exp(-h * ({inhibitor_node} - 0.5))) / ((1 - exp(0.5 * h)) * (1 + exp(-h * ({inhibitor_node} - 0.5))))) - gamma * {target_node}"
                math_ast = libsbml.parseL3Formula(formula)
                kinetic_law.setMath(math_ast)

    # Example to add plot annotations
    groups = {
        'Pro-Inflammatory': [0, 1, 2],
        'Anti-Inflammatory': [3],
        'Growth factors': [9],
        'ECM Destruction': [0, 1],
        'ECM Synthesis': [2],
        'Hypertrophy': [49]
    }
    add_plot_annotations(model, groups, node_names)

    # Write SBML file
    libsbml.writeSBMLToFile(sbml_document, filename)


filename = 'CRT.xlsx'
mact, minh, node_names, num_of_nodes = create_matrices(filename)
xinit = np.random.rand(num_of_nodes)
clamped = np.zeros(num_of_nodes)  # Example: No nodes are clamped
gamma = 1.0  # Example decay rate for all nodes
h = 10.0  # Example steepness for all nodes
export_to_sbml_with_params('model.xml', node_names, mact, minh, xinit, clamped, gamma, h)
