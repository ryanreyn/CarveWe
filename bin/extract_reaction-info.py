#!/usr/bin/env python3

import os
import argparse
from lxml import etree

def extract_ensemble_map(root):
    """Builds a dictionary mapping metaid ➝ ENSEMBLE_STATE value from RDF annotations."""
    ns = {
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
        'bqbiol': 'http://biomodels.net/biology-qualifiers/',
    }

    ensemble_map = {}
    rdf_descriptions = root.xpath("//rdf:Description", namespaces=ns)
    print(f"[DEBUG] Found {len(rdf_descriptions)} rdf:Description blocks")

    for desc in rdf_descriptions:
        about = desc.get("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about")
        if not about:
            continue
        metaid = about.replace("#", "").replace("meta_", "")
        li_elements = desc.xpath(".//rdf:li/@rdf:resource", namespaces=ns)
        for resource in li_elements:
            if "ENSEMBLE_STATE" in resource:
                raw_value = resource.split("/")[-1].strip()
                value = raw_value.replace(" ", ",")
                ensemble_map[metaid] = value
                print(f"[DEBUG] Matched ENSEMBLE_STATE: {metaid} → {value}")
                break  # Stop after first valid ENSEMBLE_STATE

    return ensemble_map

def parse_reactions_with_ensemble_states(sbml_path):
    tree = etree.parse(sbml_path)
    root = tree.getroot()

    ensemble_map = extract_ensemble_map(root)

    # Ensure SBML namespace is included
    ns = root.nsmap
    if None in ns:
        ns['sbml'] = ns.pop(None)

    reactions = root.xpath(".//sbml:reaction", namespaces=ns)
    print(f"[DEBUG] Found {len(reactions)} reaction elements")

    results = []
    for rxn in reactions:
        rxn_id = rxn.get("id") or ""
        rxn_id = rxn_id.removeprefix("R_")
        metaid_raw = rxn.get("metaid") or ""
        metaid_clean = metaid_raw.replace("meta_", "")
        ensemble_state = ensemble_map.get(metaid_clean, "")

        print(f"[DEBUG] Reaction: {rxn_id} (metaid: {metaid_raw}) → ENSEMBLE_STATE: {ensemble_state or 'None'}")
        results.append((rxn_id.strip(), ensemble_state))

    return results


def main():
    parser = argparse.ArgumentParser(description="Extract reaction IDs and ENSEMBLE_STATE values from SBML ensemble XML files.")
    parser.add_argument("-x", "--xml_dir", required=True, help="Directory containing .xml SBML files.")
    parser.add_argument("-o", "--out_dir", required=True, help="Directory to write output files.")
    args = parser.parse_args()

    xml_dir = args.xml_dir
    out_dir = args.out_dir
    output_dir = os.path.join(out_dir, "ensemble_rxn-info_files")
    filelist_path = os.path.join(out_dir, "filenames.txt")

    print("Running XML extraction:\n")

    os.makedirs(output_dir, exist_ok=True)

    xml_files = sorted(f for f in os.listdir(xml_dir) if f.endswith(".xml"))
    with open(filelist_path, "w") as f:
        for fname in xml_files:
            f.write(fname + "\n")

    for fname in xml_files:
        print(f"\n[INFO] Processing {fname}")
        prefix = os.path.splitext(fname)[0]
        out_csv = os.path.join(output_dir, f"{prefix}_rxn-info.csv")

        sbml_path = os.path.join(xml_dir, fname)
        rxn_data = parse_reactions_with_ensemble_states(sbml_path)

        with open(out_csv, "w") as f:
            for rxn_id, ensemble_state in rxn_data:
                f.write(f"{rxn_id},{ensemble_state}\n")

        print(f"[INFO] Wrote {len(rxn_data)} lines to {out_csv}")

    os.remove(filelist_path)

if __name__ == "__main__":
    main()
