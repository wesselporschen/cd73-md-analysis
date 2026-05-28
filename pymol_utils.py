from pathlib import Path

def write_view_system_frame_pml(residue_selection: list[tuple[int, str]],
                               output_path: Path,
                               frame_index: int) -> None:
    

    lines = []

    lines.append(f"load {output_path.name}")
    lines.append(f"load {(output_path.with_name(f"hull_frame_{frame_index}.pdb").name)}")
    
    lines.append(f"color cyan, frame_{frame_index} and chain A")
    lines.append(f"color green, frame_{frame_index} and chain B")
    lines.append(f"color hotpink, hull_frame_{frame_index}")

    selection = " or ".join(
            f"(resi {resi} and chain {chain})"
            for resi, chain in residue_selection
            )

    lines.append(f"select residue_selection, {selection}")
    lines.append(f"show sticks, residue_selection")
    lines.append(f"color atomic, residue_selection")
    lines.append(f"set stick_color, hotpink, residue_selection and elem C")
    lines.append(f"zoom hull_frame_{frame_index}")

    output_path.with_name(f"view_system_frame{frame_index}.pml").write_text("\n".join(lines))
    print(f"Wrote PyMol view script to {output_path.with_name(f"view_system_frame{frame_index}.pml")}")

