"""Render the four frozen signed component maps in MRIcroGL.

Launch with CPython and --mricrogl EXECUTABLE --reset-preferences.
MRIcroGL starts with -r to select radiological orientation explicitly.
This resets its saved GUI preferences. Required environment variables:

  PRONIA_COMPONENT_MAP_DIR   directory containing the four NIfTI maps
  PRONIA_PANEL_OUTPUT_DIR    destination for freshly rendered PNG panels
  PRONIA_MRICROGL_TEMPLATE   anatomical template NIfTI (for example spm152)
"""

import os


def launch(argv=None):
    import argparse
    from pathlib import Path
    import subprocess
    import time
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mricrogl", required=True, type=Path)
    parser.add_argument("--reset-preferences", required=True, action="store_true",
                        help="Allow MRIcroGL -r to reset saved GUI preferences, including radiological orientation")
    parser.add_argument("--timeout", type=float, default=300,
                        help="Maximum rendering time in seconds (default: 300)")
    args = parser.parse_args(argv)
    executable = args.mricrogl.expanduser().resolve()
    if not executable.is_file():
        parser.error("--mricrogl must name the MRIcroGL executable")
    if not 0 < args.timeout < float("inf"):
        parser.error("--timeout must be finite and positive")
    required_dir("PRONIA_COMPONENT_MAP_DIR")
    required_file("PRONIA_MRICROGL_TEMPLATE")
    if not os.environ.get("PRONIA_PANEL_OUTPUT_DIR"):
        parser.error("PRONIA_PANEL_OUTPUT_DIR must be set")
    output = Path(os.environ["PRONIA_PANEL_OUTPUT_DIR"]).expanduser().resolve()
    if output.exists():
        parser.error("PRONIA_PANEL_OUTPUT_DIR must be a new directory")
    env = dict(os.environ, PRONIA_RADIOLOGICAL_STARTUP="1")
    env["PRONIA_PANEL_OUTPUT_DIR"] = str(output)
    products = [output / "MRICROGL_RENDER_COMPLETE.txt"]
    products += [output / "comp_{}_source.png".format(i) for i in range(1, 5)]
    process = subprocess.Popen([str(executable), "-r", "-s", str(Path(__file__).resolve())], env=env)
    deadline = time.monotonic() + args.timeout
    try:
        while process.poll() is None and not products[0].is_file():
            if time.monotonic() >= deadline:
                raise RuntimeError("MRIcroGL rendering timed out before completion")
            time.sleep(0.1)
        # Some GUI builds remain open after gl.quit(). The completion marker is
        # written only after saving every panel; close only the process we started.
        try:
            code = process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            code = None
        if code not in (None, 0):
            raise RuntimeError("MRIcroGL exited with status {}".format(code))
    finally:
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
    if not all(p.is_file() and p.stat().st_size for p in products):
        raise RuntimeError("MRIcroGL did not produce all four panels and the completion marker")


def required_dir(name):
    value = os.environ.get(name, "")
    if not value or not os.path.isdir(value):
        raise RuntimeError("{} must name an existing directory".format(name))
    return os.path.abspath(value)


def required_file(name):
    value = os.environ.get(name, "")
    if not value or not os.path.isfile(value):
        raise RuntimeError("{} must name an existing file".format(name))
    return os.path.abspath(value)


def render(gl):
    # -r sets FlipLR_Radiological at application startup; resetdefaults() alone
    # preserves the user's orientation preference. Use launch() in a fresh process.
    if os.environ.get("PRONIA_RADIOLOGICAL_STARTUP") != "1":
        raise RuntimeError("Start with CPython --mricrogl EXECUTABLE --reset-preferences to set radiological orientation")
    map_dir = required_dir("PRONIA_COMPONENT_MAP_DIR")
    output_dir = os.environ.get("PRONIA_PANEL_OUTPUT_DIR", "")
    if not output_dir:
        raise RuntimeError("PRONIA_PANEL_OUTPUT_DIR must be set")
    output_dir = os.path.abspath(output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    template = required_file("PRONIA_MRICROGL_TEMPLATE")

    # Map NeuroMiner component identifiers to panel numbers.
    components = [
        ("Component_comp01_CVRabsGT3_AND_SBCGT1p3.nii", 1, "COMP01"),
        ("Component_comp02_CVRabsGT3_AND_SBCGT1p3.nii", 2, "COMP02"),
        ("Component_comp04_CVRabsGT3_AND_SBCGT1p3.nii", 3, "COMP04"),
        ("Component_comp06_CVRabsGT3_AND_SBCGT1p3.nii", 4, "COMP06"),
    ]

    # Fixed axial slice positions (mm), overlay order and signed display thresholds.
    mosaic = "A H 0.03 -50 -35 -20; -5 10 25; 40 55 70"
    positive_threshold, positive_maximum = 3.0, 15.0
    negative_threshold, negative_maximum = -3.0, -15.0

    for filename, display_number, saved_identifier in components:
        overlay = os.path.join(map_dir, filename)
        if not os.path.isfile(overlay):
            raise RuntimeError("Required signed component map is missing: " + overlay)

        gl.resetdefaults()
        gl.backcolor(255, 255, 255)
        gl.loadimage(template)
        gl.opacity(0, 60)

        gl.overlayload(overlay)
        gl.minmax(1, positive_threshold, positive_maximum)
        gl.opacity(1, 100)
        gl.colorname(1, "red")

        gl.overlayload(overlay)
        gl.minmax(2, negative_threshold, negative_maximum)
        gl.opacity(2, 80)
        gl.colorname(2, "3blue")

        gl.colorbarposition(0)
        gl.shadername("default")
        gl.linewidth(2)
        gl.bmpzoom(2)
        gl.mosaic(mosaic)

        output = os.path.join(output_dir, "comp_{}_source.png".format(display_number))
        gl.savebmp(output)
        print("Rendered Comp {} ({}) to {}".format(display_number, saved_identifier, output))

    with open(os.path.join(output_dir, "MRICROGL_RENDER_COMPLETE.txt"), "w") as marker:
        marker.write("Four component panels rendered from signed NIfTI maps.\n")
    gl.quit()


try:
    import gl
except ImportError:
    if __name__ == "__main__":
        launch()
else:
    # MRIcroGL executes scripts in its embedded namespace, which need not be
    # named __main__. Ordinary CPython imports remain inert.
    render(gl)
