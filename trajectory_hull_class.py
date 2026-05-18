from pathlib import Path
from dataclasses import dataclass
from functools import cached_property

import MDAnalysis as mda
from MDAnalysis.analysis import rms
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

from utils import write_view_system_frame_pml

@dataclass
class ConvexHullFrame:
    frame: int
    volume: float
    coords: np.ndarray
    vertex_indices: np.ndarray
    simplices: np.ndarray

class ConvexHullTrajectory():
    """
    path_to_trajectories: Path to directory containing trajectory files (e.g. .xtc), must be compatible with the topology file
    path_to_topology: Path to topology file (e.g. .pdb, .gro) compatible with the trajectory files
    residue_selection: List of tuples specifying residues to include in convex hull calculation, where each tuple is (residue_id, segment_id)
    """
    def __init__(self, path_to_trajectories: Path, path_to_topology: Path, residue_selection: list[tuple[int, str]]):
        self.path_to_trajectories = path_to_trajectories
        self.path_to_topology = path_to_topology
        self.residue_selection = residue_selection
        self._trajectories = sorted([x for x in path_to_trajectories.iterdir() if x.match('*.xtc')])
        self.universe = mda.Universe(self.path_to_topology, self._trajectories)
        self.residue_atoms = self.universe.select_atoms(" or ".join([f"(resid {resid} and segid {segid})" for resid, segid in self.residue_selection])).select_atoms("not name H*")
        self.residue_atoms_only_sidechains = self.residue_atoms.select_atoms("not backbone")
    
    @cached_property
    def hulls(self) -> list[ConvexHullFrame]:
        """
        Compute the convex hull formed from selected side chain heavy atoms for each frame of the trajectory and store the results in a list.
        """
        hulls = []

        for ts in tqdm(self.universe.trajectory):
            coords = self.residue_atoms.positions
            hull = ConvexHull(coords)
            hulls.append(
                 ConvexHullFrame(
                            frame = ts,
                            volume = hull.volume,
                            coords = coords,
                            vertex_indices = hull.vertices.copy(),
                            simplices = hull.simplices.copy()
                 )
            )

        return hulls
    
    @cached_property
    def rmsd(self) -> np.ndarray:
        """
        Compute RMSD over time for the selected residues, including all heavy atoms (side chains + backbone).
        """
        R = rms.RMSD(self.residue_atoms,
                     self.residue_atoms,
                     ref_frame=0)
        
        R.run()
        return R.rmsd[:,2]
    
    @cached_property
    def rmsd_only_sidechains(self) -> np.ndarray:
        """
        Compute RMSD over time for the selected residues, including only side chain atoms.
        """
        R = rms.RMSD(self.residue_atoms_only_sidechains,
                     self.residue_atoms_only_sidechains,
                     ref_frame=0)
        
        R.run()
        return R.rmsd[:,2]

    def plot_volume_over_time(self, title: str) -> None:
        # Raw pocket volumes
        hull_volumes = [hull.volume for hull in self.hulls]

        # Running average (1 ns = 100 frames at 10 ps/frame)
        window = 100
        average_1ns = [
            np.mean(hull_volumes[i:i+window])
            for i in range(len(hull_volumes) - window + 1)
        ]

        # Frame indices
        frames_raw = np.arange(len(hull_volumes))
        frames_avg = np.arange(window - 1, len(hull_volumes))
        frame_to_ns = lambda frame: frame * 10 / 1000 # 10 ps per frame, convert to ns
        raw_frames_as_ns = np.array([frame_to_ns(i) for i in frames_raw])
        avg_frames_as_ns = np.array([frame_to_ns(i) for i in frames_avg])


        # Seaborn styling
        sns.set_theme(style="white", context="talk")

        # Plot
        plt.figure(figsize=(12, 6))

        # Raw data: transparent
        sns.lineplot(
            x=raw_frames_as_ns,
            y=hull_volumes,
            alpha=0.25,
            linewidth=1.2,
            label="Hull volume"
        )

        # Running average: solid
        sns.lineplot(
            x=avg_frames_as_ns,
            y=average_1ns,
            linewidth=2.5,
            label="1 ns running average"
        )

        plt.xlabel("Time (ns)")
        plt.ylabel("Volume (Å³)")
        plt.title(f"{title}")

        plt.legend()
        plt.tight_layout()
  
        plt.savefig(self.path_to_trajectories / f'{title.replace(" ", "_")}.png')
        print("Plot saved to:", self.path_to_trajectories / f'{title.replace(" ", "_")}.png')
        plt.show()
    
    def plot_rmsd_over_time(self, title: str, only_sidechains: bool = False) -> None:

        if only_sidechains:
            rmsd = self.rmsd_only_sidechains
            raw_label = "RMSD (side chains only)"
        else:
            rmsd = self.rmsd
            raw_label = "RMSD (side chains + backbone)"
        
        # Running average (1 ns = 100 frames at 10 ps/frame)
        window = 100
        average_1ns = [
            np.mean(rmsd[i:i+window])
            for i in range(len(rmsd) - window + 1)
        ]

        # Frame indices
        frames_raw = np.arange(len(rmsd))
        frames_avg = np.arange(window - 1, len(rmsd))
        frame_to_ns = lambda frame: frame * 10 / 1000 # 10 ps per frame, convert to ns
        raw_frames_as_ns = np.array([frame_to_ns(i) for i in frames_raw])
        avg_frames_as_ns = np.array([frame_to_ns(i) for i in frames_avg])


        # Seaborn styling
        sns.set_theme(style="white", context="talk")

        # Plot
        plt.figure(figsize=(12, 6))

        # Raw data: transparent
        sns.lineplot(
            x=raw_frames_as_ns,
            y=rmsd,
            alpha=0.25,
            linewidth=1.2,
            label=raw_label
        )

        # Running average: solid
        sns.lineplot(
            x=avg_frames_as_ns,
            y=average_1ns,
            linewidth=2.5,
            label="1 ns running average",
            color='seagreen'
        )

        plt.xlabel("Time (ns)")
        plt.ylabel("RMSD (Å)")
        plt.title(f"{title}")

        plt.legend()
        plt.tight_layout()

        plt.savefig(self.path_to_trajectories / f'{title.replace(" ", "_")}.png')
        print("Plot saved to:", self.path_to_trajectories / f'{title.replace(" ", "_")}.png')
        plt.show()


    def write_convex_hull_frame_pdb(self, hull_frame: ConvexHullFrame, output_path: Path | None = None) -> None:
        """
        Write a single convex hull frame to a PDB file, including vertex coordinates and connectivity.
        """
        if output_path is None:
            output_path = self.path_to_trajectories / "hull_frame.pdb"

        coords = hull_frame.coords
        vertex_indices = hull_frame.vertex_indices
        simplices = hull_frame.simplices

        vertex_set = set(vertex_indices)

        # original coord index -> pdb atom index
        index_map = {}

        with open(output_path, "w") as f:

            pdb_index = 1

            # write hull vertices only
            for original_idx in vertex_indices:

                point = coords[original_idx]

                f.write(
                    f"HETATM{pdb_index:5d}  C   BOX A   1"
                    f"{point[0]:12.3f}"
                    f"{point[1]:8.3f}"
                    f"{point[2]:8.3f}"
                    f"  1.00  0.00           C\n"
                )

                index_map[original_idx] = pdb_index
                pdb_index += 1

            # build unique edge list from simplices
            edges = set()

            for tri in simplices:

                for i, j in [(0, 1), (1, 2), (2, 0)]:

                    a = tri[i]
                    b = tri[j]

                    if a in vertex_set and b in vertex_set:

                        edge = tuple(sorted((a, b)))
                        edges.add(edge)

            # write connectivity
            for a, b in edges:

                f.write(
                    f"CONECT"
                    f"{index_map[a]:5d}"
                    f"{index_map[b]:5d}\n"
                )

            f.write("END\n")

        print(f"Wrote hull PDB to: {output_path}")

    def write_system_frame_pdb(self, frame_index: int, output_path: Path | None = None) -> None:
        """
        Write the full system for a given frame index to a PDB file and also write the corresponding convex hull frame to a separate PDB file.
        """
        if output_path is None:
            output_path = self.path_to_trajectories / f"frame_{frame_index}.pdb"

        with mda.Writer(output_path, multiframe=False) as writer:
            self.universe.trajectory[frame_index]
            writer.write(self.universe)

        print(f"Wrote system frame {frame_index} to PDB: {output_path}")

        self.write_convex_hull_frame_pdb(hull_frame=self.hulls[frame_index], 
                                         output_path=output_path.with_name(f"hull_frame_{frame_index}.pdb"))
        write_view_system_frame_pml(residue_selection=self.residue_selection, output_path=output_path, frame_index=frame_index)

