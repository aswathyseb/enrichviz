# Prepend this directory so `enrichviz` is on PATH inside pixi shell.
if [ -n "${BASH_SOURCE[0]:-}" ]; then
  _enrichviz_bin=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
else
  _enrichviz_bin="${PIXI_PROJECT_ROOT}/bin"
fi

export PATH="${_enrichviz_bin}:${PATH}"

unset _enrichviz_bin
