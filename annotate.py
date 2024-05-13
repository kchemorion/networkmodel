import pandas as pd
from lxml import etree as ET

def read_excel_annotations(file_path):
    # Namespace map for RDF and other namespaces
    ns_map = {
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
        'dc': 'http://purl.org/dc/elements/1.1/',
        'vCard': 'http://www.w3.org/2001/vcard-rdf/3.0#',
    }
    
    # Read data from Excel sheets
    references = pd.read_excel(file_path, sheet_name='References')
    creators = pd.read_excel(file_path, sheet_name='Creators')
    plot_config_data = pd.read_excel(file_path, sheet_name='PlotConfig')
    print("Columns in PlotConfig:", plot_config_data.columns) 
    
    # Create XML elements for RDF metadata
    rdf_root = ET.Element('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}RDF', nsmap=ns_map)
    
    # Add creator details to RDF
    desc = ET.SubElement(rdf_root, '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}Description', attrib={'{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about': '#ModelCreators'})
    creator_bag = ET.SubElement(desc, '{http://purl.org/dc/elements/1.1/}creator', attrib={'{http://www.w3.org/1999/02/22-rdf-syntax-ns#}parseType': "Resource"})
    
    for _, creator in creators.iterrows():
        li = ET.SubElement(creator_bag, '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}li', attrib={'{http://www.w3.org/1999/02/22-rdf-syntax-ns#}parseType': "Resource"})
        n = ET.SubElement(li, '{http://www.w3.org/2001/vcard-rdf/3.0#}N', attrib={'{http://www.w3.org/1999/02/22-rdf-syntax-ns#}parseType': "Resource"})
        ET.SubElement(n, '{http://www.w3.org/2001/vcard-rdf/3.0#}Family').text = creator['FamilyName']
        ET.SubElement(n, '{http://www.w3.org/2001/vcard-rdf/3.0#}Given').text = creator['GivenName']
        ET.SubElement(li, '{http://www.w3.org/2001/vcard-rdf/3.0#}EMAIL').text = creator['Email']
        org = ET.SubElement(li, '{http://www.w3.org/2001/vcard-rdf/3.0#}ORG', attrib={'{http://www.w3.org/1999/02/22-rdf-syntax-ns#}parseType': "Resource"})
        ET.SubElement(org, '{http://www.w3.org/2001/vcard-rdf/3.0#}Orgname').text = creator['Organization']

    # Handle references
    desc_references = ET.SubElement(rdf_root, '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}Description', attrib={'{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about': '#ModelReferences'})
    for _, ref in references.iterrows():
        for i in range(1, 4):  # Assuming you have 3 reference columns
            reference_col = f'Reference{i}'
            if pd.notna(ref[reference_col]):
                ET.SubElement(desc_references, '{http://purl.org/dc/elements/1.1/}source').text = ref[reference_col]
    
    # Create plot configuration annotations
    plot_config = ET.Element('plot_config')
    
    for _, row in plot_config_data.iterrows():
        group_elem = ET.SubElement(plot_config, 'plot_group', attrib={'name': row['Group'], 'type': row['Type']})
        species_ids = row['SpeciesID'].split(',')
        for species_id in species_ids:
            ET.SubElement(group_elem, 'species', attrib={'id': species_id.strip()})
    
    # Convert XML to string
    creators_xml_string = ET.tostring(rdf_root, pretty_print=True, xml_declaration=True, encoding='UTF-8')
    plot_config_xml_string = ET.tostring(plot_config, pretty_print=True, xml_declaration=True, encoding='UTF-8')

    # Print XML strings
    print(creators_xml_string.decode('utf-8'))
    print(plot_config_xml_string.decode('utf-8'))
    
# Example usage
    xml_annotation = read_excel_annotations('annotations.xlsx')
