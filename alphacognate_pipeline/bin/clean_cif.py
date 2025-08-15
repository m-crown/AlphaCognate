import sys


def filter_cif_by_plddt(input_filepath, output_filepath, plddt_min_threshold=70.0):
    """
    Filters a CIF file to remove residues with pLDDT below a threshold.
    """
    try:
        with open(input_filepath, 'r') as f_in:
            lines = f_in.readlines()
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_filepath}")
        return
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    atom_site_loop_start_idx = -1
    atom_site_header_lines = []

    for i, line in enumerate(lines):
        if line.strip() == "loop_":
            if i + 1 < len(lines) and lines[i + 1].strip().startswith("_atom_site."):
                atom_site_loop_start_idx = i
                break

    if atom_site_loop_start_idx == -1:
        print("Warning: '_atom_site' loop not found. Writing original content to output.")
        try:
            with open(output_filepath, 'w') as f_out:
                f_out.writelines(lines)
        except Exception as e:
            print(f"Error writing output file: {e}")
        return

    header_parsing_idx = atom_site_loop_start_idx + 1
    while header_parsing_idx < len(lines) and lines[header_parsing_idx].strip().startswith("_atom_site."):
        atom_site_header_lines.append(lines[header_parsing_idx])
        header_parsing_idx += 1

    atom_data_actual_start_idx = header_parsing_idx

    atom_data_lines_raw_for_parsing = []
    atom_site_data_block_end_idx = atom_data_actual_start_idx

    for i in range(atom_data_actual_start_idx, len(lines)):
        line = lines[i]
        line_stripped = line.strip()

        if not line_stripped:
            atom_site_data_block_end_idx = i
            break

        if line_stripped.startswith(("_", "loop_", "data_", "save_", "#")):
            atom_site_data_block_end_idx = i
            break

        atom_data_lines_raw_for_parsing.append(line)
        atom_site_data_block_end_idx = i + 1

    if not atom_site_header_lines or not atom_data_lines_raw_for_parsing:
        print("Warning: Could not parse _atom_site headers or data. Writing original content.")
        try:
            with open(output_filepath, 'w') as f_out:
                f_out.writelines(lines)
        except Exception as e:
            print(f"Error writing output file: {e}")
        return

    column_names = [h.strip() for h in atom_site_header_lines]
    try:
        plddt_col_idx = column_names.index("_atom_site.B_iso_or_equiv")
        chain_id_col_idx = column_names.index("_atom_site.label_asym_id")
        res_seq_col_idx = column_names.index("_atom_site.label_seq_id")
    except ValueError:
        print(
            "Error: Required columns (_atom_site.B_iso_or_equiv, _atom_site.label_asym_id, _atom_site.label_seq_id) not all found in _atom_site loop.")
        return

    residues_to_remove = set()
    parsed_atom_data_with_info = []

    for data_line_raw in atom_data_lines_raw_for_parsing:
        fields = data_line_raw.strip().split()

        if len(fields) != len(column_names):
            parsed_atom_data_with_info.append({
                "line": data_line_raw, "chain": None, "res_seq": None, "valid_plddt_info": False
            })
            continue

        try:
            plddt_str = fields[plddt_col_idx]
            chain_id = fields[chain_id_col_idx]
            res_seq_num = fields[res_seq_col_idx]

            current_atom_info = {
                "line": data_line_raw,
                "chain": chain_id,
                "res_seq": res_seq_num,
                "valid_plddt_info": True
            }
            parsed_atom_data_with_info.append(current_atom_info)

            plddt = float(plddt_str)
            if plddt < plddt_min_threshold:
                residues_to_remove.add((chain_id, res_seq_num))

        except (IndexError, ValueError):
            parsed_atom_data_with_info[-1]["valid_plddt_info"] = False

    try:
        with open(output_filepath, 'w') as f_out:
            f_out.writelines(lines[:atom_site_loop_start_idx])
            f_out.write(lines[atom_site_loop_start_idx])
            for header_line in atom_site_header_lines:
                f_out.write(header_line)

            for atom_entry in parsed_atom_data_with_info:
                if atom_entry["valid_plddt_info"]:
                    if (atom_entry["chain"], atom_entry["res_seq"]) not in residues_to_remove:
                        f_out.write(atom_entry["line"])
                else:
                    f_out.write(atom_entry["line"])
            f_out.writelines(lines[atom_site_data_block_end_idx:])
        print(f"Filtered CIF file written to {output_filepath}")
    except Exception as e:
        print(f"Error writing output file: {e}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_cif_file> <output_cif_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    filter_cif_by_plddt(input_file, output_file)

    # Example usage within Python:
    # filter_cif_by_plddt("input.cif", "output_filtered.cif", plddt_min_threshold=70.0)