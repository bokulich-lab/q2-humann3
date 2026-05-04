# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pathlib import Path
import shutil
import subprocess


EXTERNAL_CMD_WARNING = (
    "Running external command line application. This may print additional "
    "output below."
)


def run_command(
    cmd: list[str], env=None, verbose: bool = True, pipe: bool = False,
    **kwargs
):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(" ".join(cmd), end="\n\n")

    if pipe:
        result = subprocess.run(
            cmd,
            env=env,
            check=True,
            capture_output=True,
            text=True,
            **kwargs,
        )
        return result

    if env:
        subprocess.run(cmd, env=env, check=True, **kwargs)
    else:
        subprocess.run(cmd, check=True, **kwargs)


def run_humann_command(cmd: list[str]) -> None:
    try:
        run_command(cmd)
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        stdout = (exc.stdout or "").strip()
        detail = stderr or stdout or str(exc)
        raise RuntimeError(
            f"Command failed with exit code {exc.returncode}: {detail}"
        ) from exc


def copy_directory_contents(source_dir: Path, target_dir: Path) -> None:
    for path in source_dir.iterdir():
        destination = target_dir / path.name
        if path.is_dir():
            shutil.copytree(path, destination)
        else:
            shutil.copy2(path, destination)
