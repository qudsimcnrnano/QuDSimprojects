

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import subprocess


class LoadFileHandler:
    """Manages file loading operations for the QuDSim GUI."""

    FILETYPES = {
        "geo": [("GMSH Geometry", "*.geo"), ("All files", "*.*")],
        "msh": [("GMSH Mesh", "*.msh"), ("All files", "*.*")],
        "vtu": [("VTK Unstructured Grid", "*.vtu"),
                ("VTK files", "*.vtk"), ("All files", "*.*")],
    }

    def __init__(self, app):
        self.app = app

    def load_file(self, file_type="geo"):
        """
        Open a file dialog and load the selected file.

        Parameters
        ----------
        file_type : str
            One of 'geo', 'msh', or 'vtu'.
        """
        filetypes = self.FILETYPES.get(file_type, self.FILETYPES["geo"])
        filepath = filedialog.askopenfilename(
            title="Open .{} File".format(file_type),
            filetypes=filetypes)

        if not filepath:
            return

        self.app.update_status("Loading {}...".format(
            os.path.basename(filepath)))

        if file_type == "geo":
            self._load_geo(filepath)
        elif file_type == "msh":
            self._load_msh(filepath)
        elif file_type == "vtu":
            self._load_vtu(filepath)

    # ------------------------------------------------------------------ #
    #  .geo files — open in text editor                                    #
    # ------------------------------------------------------------------ #
    def _load_geo(self, filepath):
        """Load and display a .geo file in an editor window."""
        self.app.current_geo_file = filepath

        editor = tk.Toplevel(self.app.root)
        editor.title("Geometry Editor - {}".format(os.path.basename(filepath)))
        editor.geometry("750x550")

        # Toolbar
        toolbar = tk.Frame(editor, bg="#eef2f7", padx=5, pady=3)
        toolbar.pack(fill=tk.X)

        tk.Button(toolbar, text="Save",
                  command=lambda: self._save_geo(filepath, text_area),
                  font=("Helvetica", 9, "bold")).pack(side=tk.LEFT, padx=3)
        tk.Button(toolbar, text="Generate Mesh",
                  command=self.app._run_gmsh,
                  font=("Helvetica", 9, "bold"),
                  bg="#2ecc71", fg="white").pack(side=tk.LEFT, padx=3)
        tk.Label(toolbar, text=filepath, font=("Helvetica", 8),
                 fg="#5d6d7e", bg="#eef2f7").pack(side=tk.RIGHT, padx=5)

        # Text editor with syntax highlighting colours
        text_frame = tk.Frame(editor)
        text_frame.pack(fill=tk.BOTH, expand=True)

        # Line numbers
        line_numbers = tk.Text(text_frame, width=5, padx=4,
                               bg="#2c3e50", fg="#95a5a6",
                               font=("Courier", 11),
                               state=tk.DISABLED, takefocus=0)
        line_numbers.pack(side=tk.LEFT, fill=tk.Y)

        text_area = tk.Text(text_frame, font=("Courier", 11),
                            wrap=tk.NONE, undo=True, bg="#fafafa",
                            insertbackground="#2c3e50")
        ys = ttk.Scrollbar(text_frame, orient=tk.VERTICAL,
                           command=text_area.yview)
        xs = ttk.Scrollbar(editor, orient=tk.HORIZONTAL,
                           command=text_area.xview)
        text_area.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)

        ys.pack(side=tk.RIGHT, fill=tk.Y)
        text_area.pack(fill=tk.BOTH, expand=True)
        xs.pack(fill=tk.X)

        # Load content
        try:
            with open(filepath, "r") as f:
                content = f.read()
            text_area.insert(tk.END, content)
            self.app.update_status("Loaded: {}".format(
                os.path.basename(filepath)))
        except Exception as e:
            messagebox.showerror("Load Error",
                                 "Failed to load file:\n{}".format(str(e)))

        # Simple syntax colouring for GMSH keywords
        self._apply_geo_highlighting(text_area)

        # Update line numbers
        def update_lines(*args):
            line_numbers.config(state=tk.NORMAL)
            line_numbers.delete("1.0", tk.END)
            line_count = int(text_area.index("end-1c").split(".")[0])
            line_numbers.insert(tk.END,
                                "\n".join(str(i) for i in range(1, line_count + 1)))
            line_numbers.config(state=tk.DISABLED)

        text_area.bind("<KeyRelease>", update_lines)
        update_lines()

    def _save_geo(self, filepath, text_widget):
        content = text_widget.get("1.0", tk.END)
        with open(filepath, "w") as f:
            f.write(content)
        self.app.update_status("Saved: {}".format(os.path.basename(filepath)))

    def _apply_geo_highlighting(self, text_widget):
        """Basic keyword highlighting for .geo files."""
        keywords = [
            "Point", "Line", "Circle", "Plane Surface", "Surface",
            "Volume", "Physical", "Extrude", "Mesh", "Transfinite",
            "Recombine", "SetFactory", "BooleanDifference",
            "BooleanUnion", "BooleanIntersection",
        ]
        text_widget.tag_configure("keyword", foreground="#2e86c1",
                                  font=("Courier", 11, "bold"))
        text_widget.tag_configure("comment", foreground="#27ae60")
        text_widget.tag_configure("number", foreground="#e67e22")

        content = text_widget.get("1.0", tk.END)
        for i, line in enumerate(content.split("\n"), start=1):
            # Comments
            if "//" in line:
                idx = line.index("//")
                text_widget.tag_add("comment",
                                    "{}.{}".format(i, idx),
                                    "{}.end".format(i))
            # Keywords
            for kw in keywords:
                start = 0
                while True:
                    pos = line.find(kw, start)
                    if pos == -1:
                        break
                    text_widget.tag_add("keyword",
                                        "{}.{}".format(i, pos),
                                        "{}.{}".format(i, pos + len(kw)))
                    start = pos + len(kw)

    # ------------------------------------------------------------------ #
    #  .msh files — visualize in GMSH or show info                         #
    # ------------------------------------------------------------------ #
    def _load_msh(self, filepath):
        """Load a .msh file and offer visualization options."""
        self.app.current_msh_file = filepath

        # Try to get mesh statistics
        info = self._get_mesh_info(filepath)

        # Display mesh info window
        info_win = tk.Toplevel(self.app.root)
        info_win.title("Mesh Info - {}".format(os.path.basename(filepath)))
        info_win.geometry("500x400")

        tk.Label(info_win, text="Mesh File Information",
                 font=("Helvetica", 14, "bold"),
                 fg="#2c3e50").pack(pady=(15, 5))
        tk.Label(info_win, text=os.path.basename(filepath),
                 font=("Helvetica", 10),
                 fg="#7f8c8d").pack(pady=(0, 10))

        # Info display
        info_text = tk.Text(info_win, font=("Courier", 11),
                            bg="#fafafa", height=12, wrap=tk.WORD)
        info_text.pack(fill=tk.BOTH, expand=True, padx=15, pady=5)
        info_text.insert(tk.END, info)
        info_text.config(state=tk.DISABLED)

        # Action buttons
        btn_frame = tk.Frame(info_win, pady=10)
        btn_frame.pack(fill=tk.X, padx=15)

        tk.Button(btn_frame, text="Open in GMSH",
                  font=("Helvetica", 10, "bold"),
                  bg="#2ecc71", fg="white",
                  command=lambda: self._open_in_gmsh(filepath)).pack(
                      side=tk.LEFT, padx=5)
        tk.Button(btn_frame, text="Close",
                  font=("Helvetica", 10),
                  command=info_win.destroy).pack(side=tk.RIGHT, padx=5)

        self.app.update_status("Loaded mesh: {}".format(
            os.path.basename(filepath)))

    def _get_mesh_info(self, filepath):
        """Parse basic mesh statistics from a .msh file."""
        nodes = 0
        elements = 0
        dim = "Unknown"
        file_size = os.path.getsize(filepath) / 1024  # KB

        try:
            with open(filepath, "r") as f:
                section = None
                for line in f:
                    line = line.strip()
                    if line == "$Nodes":
                        section = "nodes"
                        continue
                    elif line == "$EndNodes":
                        section = None
                        continue
                    elif line == "$Elements":
                        section = "elements"
                        continue
                    elif line == "$EndElements":
                        section = None
                        continue

                    if section == "nodes" and nodes == 0:
                        try:
                            nodes = int(line.split()[0])
                        except (ValueError, IndexError):
                            pass
                    elif section == "elements" and elements == 0:
                        try:
                            elements = int(line.split()[0])
                        except (ValueError, IndexError):
                            pass
        except Exception:
            pass

        info = (
            "File:       {}\n"
            "Size:       {:.1f} KB\n"
            "Nodes:      {}\n"
            "Elements:   {}\n"
            "\n"
            "This mesh file can be visualized in GMSH\n"
            "or used as input for QuDSim simulations.\n"
        ).format(os.path.basename(filepath), file_size, nodes, elements)

        return info

    def _open_in_gmsh(self, filepath):
        """Open mesh file in GMSH for visualization."""
        try:
            subprocess.Popen(["gmsh", filepath])
            self.app.update_status("Opened in GMSH")
        except FileNotFoundError:
            messagebox.showerror(
                "GMSH Not Found",
                "GMSH is not installed or not in PATH.\n"
                "Install from: https://gmsh.info")

    # ------------------------------------------------------------------ #
    #  .vtu files — open in ParaView                                       #
    # ------------------------------------------------------------------ #
    def _load_vtu(self, filepath):
        """Load a .vtu file and open in ParaView."""
        self.app.current_vtu_file = filepath
        self.app.update_status("Opening {} in ParaView...".format(
            os.path.basename(filepath)))
        try:
            subprocess.Popen(["paraview", filepath])
            self.app.update_status("ParaView launched")
        except FileNotFoundError:
            # Fall back to pvpython if available
            try:
                subprocess.Popen(["pvpython", "--", filepath])
                self.app.update_status("pvpython launched")
            except FileNotFoundError:
                messagebox.showerror(
                    "ParaView Not Found",
                    "Neither paraview nor pvpython found in PATH.\n"
                    "Install ParaView from: https://www.paraview.org")
