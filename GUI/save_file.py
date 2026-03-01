
import tkinter as tk
from tkinter import filedialog, messagebox
import os
import json
from datetime import datetime


class SaveFileHandler:
    """Manages file saving operations for the QuDSim GUI."""

    def __init__(self, app):
        self.app = app

    def save_current(self):
        """Save the current project configuration."""
        config = {
            "project_name":     self.app.project_name.get(),
            "modules_required": self.app.modules_required.get(),
            "project_version":  self.app.project_version.get(),
            "maintainer_email": self.app.maintainer_email.get(),
            "mpi_processes":    self.app.mpi_processes.get(),
            "geo_file":         self.app.current_geo_file,
            "msh_file":         self.app.current_msh_file,
            "vtu_file":         self.app.current_vtu_file,
            "input_param_file": self.app.input_param_file,
            "saved_at":         datetime.now().isoformat(),
        }

        # Determine save path
        project_name = self.app.project_name.get().strip()
        if project_name:
            default_name = project_name.replace(" ", "_") + ".qdsproj"
        else:
            default_name = "untitled.qdsproj"

        filepath = os.path.join(os.getcwd(), default_name)

        if not os.path.exists(filepath):
            filepath = self.save_as()
            if not filepath:
                return
        else:
            try:
                with open(filepath, "w") as f:
                    json.dump(config, f, indent=2)
                self.app.update_status(
                    "Project saved: {}".format(os.path.basename(filepath)))
            except Exception as e:
                messagebox.showerror("Save Error", str(e))

    def save_as(self):
        """Save project configuration to a user-selected location."""
        project_name = self.app.project_name.get().strip()
        if project_name:
            default_name = project_name.replace(" ", "_") + ".qdsproj"
        else:
            default_name = "untitled.qdsproj"

        filepath = filedialog.asksaveasfilename(
            title="Save Project As",
            initialfile=default_name,
            defaultextension=".qdsproj",
            filetypes=[("QuDSim Project", "*.qdsproj"),
                       ("JSON files", "*.json"),
                       ("All files", "*.*")])

        if not filepath:
            return None

        config = {
            "project_name":     self.app.project_name.get(),
            "modules_required": self.app.modules_required.get(),
            "project_version":  self.app.project_version.get(),
            "maintainer_email": self.app.maintainer_email.get(),
            "mpi_processes":    self.app.mpi_processes.get(),
            "geo_file":         self.app.current_geo_file,
            "msh_file":         self.app.current_msh_file,
            "vtu_file":         self.app.current_vtu_file,
            "input_param_file": self.app.input_param_file,
            "saved_at":         datetime.now().isoformat(),
        }

        try:
            with open(filepath, "w") as f:
                json.dump(config, f, indent=2)
            self.app.update_status(
                "Project saved: {}".format(os.path.basename(filepath)))
            return filepath
        except Exception as e:
            messagebox.showerror("Save Error", str(e))
            return None

    def export_geo(self, text_widget, filepath=None):
        """Export the contents of a text widget as a .geo file."""
        if not filepath:
            filepath = filedialog.asksaveasfilename(
                title="Export .geo File",
                defaultextension=".geo",
                filetypes=[("GMSH Geometry", "*.geo"),
                           ("All files", "*.*")])
        if not filepath:
            return

        content = text_widget.get("1.0", tk.END)
        try:
            with open(filepath, "w") as f:
                f.write(content)
            self.app.current_geo_file = filepath
            self.app.update_status(
                "Exported: {}".format(os.path.basename(filepath)))
        except Exception as e:
            messagebox.showerror("Export Error", str(e))
