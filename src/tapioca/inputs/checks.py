'''
Checks of the input parameters of the scenario builders (tapioca.inputs).

InputChecks adds standardized checks to a class, as VerboseExtension adds the verbose prints: each check returns
the value checked (normalized) or raises `self.error` with a message that names the parameter, so the classes
that build the scenarios do not repeat the same `if ... raise` blocks. SVGChecks extends it with the checks of
the scenarios drawn in SVG files (tapioca.inputs.svg.SVGScenario).
'''

import os

__all__ = ["InputError", "SVGInputError", "InputChecks", "SVGChecks"]


class InputError(ValueError):
    '''Error in the input parameters of a scenario'''


class SVGInputError(InputError):
    '''Error in the SVG file or in the settings of the scenario'''


class InputChecks:
    '''Adds standardized checks of the input parameters. `error` is the exception raised by the checks.'''
    error = InputError

    def _fail(self, message: str):
        raise self.error(message)

    def _warn(self, message: str):
        print(f"Warning: {message}")

    def _check_choice(self, name: str, value, choices: list, aliases: dict = None) -> str:
        '''
        Check that value is one of choices (case insensitive) or of the aliases ({alias: choice}).
        Return the choice as written in choices.
        '''
        key = str(value).strip()
        options = {choice.lower(): choice for choice in choices}
        options.update({alias.lower(): choice for alias, choice in (aliases or {}).items()})
        if key.lower() not in options:
            self._fail(f"Unknown {name} '{key}'. Use {', '.join(choices)}.")
        return options[key.lower()]

    def _check_on_off(self, name: str, value) -> str:
        '''Check that value is 'on' or 'off'.'''
        return self._check_choice(name, value, ['on', 'off'])

    def _check_number(self, name: str, value, minimum: float = None, strict: bool = False) -> float:
        '''Check that value is a number, at least minimum (greater than minimum if strict).'''
        try:
            value = float(value)
        except (TypeError, ValueError):
            self._fail(f"{name} must be a number (found '{value}').")
        if minimum is not None and (value < minimum or (strict and value == minimum)):
            self._fail(f"{name} must be {'greater than' if strict else 'at least'} {minimum:g} (found {value:g}).")
        return value

    def _check_path(self, name: str, path, base: str = None) -> str:
        '''Check that the file or folder path exists (a relative path is relative to the folder base, if given).'''
        path = os.path.expanduser(str(path).strip())
        if base is not None and not os.path.isabs(path):
            path = os.path.join(base, path)
        path = os.path.normpath(path)
        if not os.path.exists(path):
            self._fail(f"{name} not found: {path}")
        return path

    def _check_file(self, name: str, path) -> str:
        '''Check that the file exists.'''
        if not os.path.isfile(path):
            self._fail(f"{name} not found: {path}")
        return path

    def _check_keys(self, name: str, keys, allowed) -> None:
        '''Check that all keys are in allowed.'''
        unknown = sorted(set(keys) - set(allowed))
        if unknown:
            self._fail(f"Unknown {name}: {', '.join(unknown)}. Use {', '.join(allowed)}.")

    def _check_domain(self, Lx: float, Lz: float, Nx: int, Nz: int) -> None:
        '''Check the size (m) and the number of nodes of a 2D domain.'''
        self._check_number('Lx', Lx, 0, strict=True)
        self._check_number('Lz', Lz, 0, strict=True)
        self._check_number('Nx', Nx, 2)
        self._check_number('Nz', Nz, 2)


class SVGChecks(InputChecks):
    '''Checks of the scenarios drawn in SVG files (tapioca.inputs.svg.SVGScenario).'''
    error = SVGInputError

    def _check_setup(self, setup: dict, svg_file: str, label: str) -> dict:
        '''Check that the SVG has the setup shape (read by tapioca.inputs.svg.read_setup).'''
        if setup is None:
            self._fail(f"No shape labelled '{label}' in {svg_file}. Draw a rectangle over the model domain, "
                       f"label it '{label}' and write Lx, Lz, dx and dz in its description.")
        return setup

    def _check_svg_versions(self, svg_file: str, svg_version: str, inkscape_version: str,
                            tested_svg: str, tested_inkscape: str) -> None:
        '''Warn when the SVG declares another version of the format or was saved by another version of Inkscape.'''
        if svg_version is not None and svg_version != tested_svg:
            self._warn(f"{svg_file} declares SVG {svg_version}; tapioca reads SVG {tested_svg}. Check the layers.")
        if inkscape_version is not None and not inkscape_version.startswith(tested_inkscape + "."):
            self._warn(f"{svg_file} was saved by Inkscape {inkscape_version.split()[0]}; the SVG scenarios were tested "
                       f"with Inkscape {tested_inkscape}. Check the layers (e.g. tapioca-svg-interfaces --plot).")

    def _check_layers(self, layers: list, svg_file: str, asthenosphere_label: str) -> None:
        '''Check the stack of layers (bottom to top): an asthenosphere and the sticky air at the top.'''
        if layers[-1].rheology_name != 'air':
            self._fail(f"The top layer of {svg_file} must be the sticky air, with 'air' in the first line of its "
                       f"description (found '{layers[-1].layer_label}' with {layers[-1].rheology_name}).")
        if not any(layer.layer_label == asthenosphere_label for layer in layers):
            self._fail(f"No shape labelled '{asthenosphere_label}' in {svg_file}.")

    def _check_experiment(self, regime: str, moving_boundary: str, source: dict, magmatism: str,
                          magmatism_extraction: str) -> None:
        '''Check the combinations of settings of the numerical experiment.'''
        if moving_boundary == 'symmetric' and regime == 'compression':
            self._fail("moving_boundary = symmetric is only available in the rifting and accordion regimes.")
        if regime in ['rifting', 'compression'] and source.get('quiescence_time') != 'default':
            self._warn(f"quiescence_time is not used in the {regime} regime.")
        if magmatism_extraction == 'on' and magmatism == 'off':
            self._warn("magmatism_extraction is on but magmatism is off.")
