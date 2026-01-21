#!/usr/bin/env python3

import os
import argparse
from lxml import etree

def extract_ensemble_map(root):
    """Builds a dictionary mapping reaction_id/metaid ➝ ENSEMBLE_STATE value.

    Supports both:
    - Legacy: Notes section HTML (mapped by reaction ID)
    - New: RDF annotations (mapped by metaid)
    """
    ensemble_map = {}

    # Method 1: Try notes section (legacy SBML format)
    ns = root.nsmap.copy()
    if None in ns:
        ns['sbml'] = ns.pop(None)
    ns['html'] = 'http://www.w3.org/1999/xhtml'

    reactions = root.xpath(".//sbml:reaction", namespaces=ns)
    notes_found = 0

    for rxn in reactions:
        rxn_id = rxn.get("id", "")
        if not rxn_id:
            continue

        rxn_id_clean = rxn_id.replace("R_", "")

        # Look for ENSEMBLE_STATE in notes
        p_tags = rxn.xpath(".//html:p", namespaces=ns) or rxn.xpath(".//p")

        for p in p_tags:
            p_text = (p.text or "").strip()
            if p_text.startswith('ENSEMBLE_STATE:'):
                state_str = p_text.replace('ENSEMBLE_STATE:', '').strip()
                state_csv = state_str.replace(' ', ',')
                ensemble_map[rxn_id_clean] = state_csv
                notes_found += 1
                print(f"[DEBUG] Notes (legacy): {rxn_id_clean} → {state_csv[:50]}...")
                break

    if notes_found > 0:
        print(f"[DEBUG] Found {notes_found} ensemble states in notes (legacy format)")

    # Method 2: Try RDF annotations (new SBML format)
    ns_rdf = {
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
        'bqbiol': 'http://biomodels.net/biology-qualifiers/',
    }

    rdf_descriptions = root.xpath("//rdf:Description", namespaces=ns_rdf)
    if rdf_descriptions:
        print(f"[DEBUG] Found {len(rdf_descriptions)} RDF blocks (new format)")
        for desc in rdf_descriptions:
            about = desc.get("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about")
            if not about:
                continue
            metaid = about.replace("#", "").replace("meta_", "")
            li_elements = desc.xpath(".//rdf:li/@rdf:resource", namespaces=ns_rdf)
            for resource in li_elements:
                if "ENSEMBLE_STATE" in resource:
                    raw_value = resource.split("/")[-1].strip()
                    value = raw_value.replace(" ", ",")
                    ensemble_map[metaid] = value
                    print(f"[DEBUG] RDF (new): {metaid} → {value[:50]}...")
                    break

    print(f"[DEBUG] Total ensemble states extracted: {len(ensemble_map)}")
    return ensemble_map

def parse_reactions_with_ensemble_states(sbml_path):
    try:
        tree = etree.parse(sbml_path)
    except etree.XMLSyntaxError as e:
        # Attempt to clean null bytes and retry
        print(f"[WARN] XML parsing error in {sbml_path}: {e}")
        print(f"[INFO] Attempting to clean null bytes and retry...")
        try:
            with open(sbml_path, 'rb') as f:
                content = f.read()
            content_clean = content.replace(b'\x00', b'')

            # Parse from cleaned bytes
            tree = etree.fromstring(content_clean)
            root = tree
        except Exception as clean_error:
            print(f"[ERROR] Failed to parse even after cleaning: {clean_error}")
            raise
    else:
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
        rxn_id_clean = rxn_id.removeprefix("R_")
        metaid_raw = rxn.get("metaid") or ""
        metaid_clean = metaid_raw.replace("meta_", "").replace("R_", "")

        # Try lookup by reaction ID first (legacy notes), then by metaid (new RDF)
        ensemble_state = ensemble_map.get(rxn_id_clean) or ensemble_map.get(metaid_clean, "")

        status = "Found" if ensemble_state else "None"
        print(f"[DEBUG] Reaction: {rxn_id_clean} → ENSEMBLE_STATE: {status}")
        results.append((rxn_id_clean.strip(), ensemble_state))

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

        try:
            rxn_data = parse_reactions_with_ensemble_states(sbml_path)

            with open(out_csv, "w") as f:
                for rxn_id, ensemble_state in rxn_data:
                    f.write(f"{rxn_id},{ensemble_state}\n")

            print(f"[INFO] Wrote {len(rxn_data)} lines to {out_csv}")
        except Exception as e:
            print(f"[ERROR] Failed to process {fname}: {e}")
            print(f"[WARN] Skipping {fname} and continuing with next file...")
            continue

    os.remove(filelist_path)

if __name__ == "__main__":
    main()
