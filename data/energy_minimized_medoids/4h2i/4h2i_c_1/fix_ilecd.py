with open("nosol.pdb") as f, open("nosol_rename.pdb", "w") as out:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20]
            atomname = line[12:16]

            if resname == "ILE" and atomname == " CD ":
                line = line[:12] + " CD1" + line[16:]

        out.write(line)
