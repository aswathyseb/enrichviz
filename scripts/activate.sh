# Prepend the repo bin directory so `enrichviz` is on PATH in pixi shell / pixi run.
if [ -n "${BASH_SOURCE[0]}" ]; then
  _enrichviz_scripts=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
else
  _enrichviz_scripts="${PIXI_PROJECT_ROOT}/scripts"
fi
_enrichviz_root=$(CDPATH= cd -- "${_enrichviz_scripts}/.." && pwd)
_enrichviz_comp="${_enrichviz_scripts}/enrichviz-completion.sh"

export PATH="${_enrichviz_root}/bin:${PATH}"

# `pixi shell` starts a new bash with exported env only. `complete` does not
# survive that, so install a PROMPT_COMMAND hook the new shell will run.
case "${PROMPT_COMMAND:-}" in
  *enrichviz-completion.sh*) ;;
  *)
    export PROMPT_COMMAND=". '${_enrichviz_comp}'${PROMPT_COMMAND:+; ${PROMPT_COMMAND}}"
    ;;
esac

# When this file is sourced in the current shell (pixi shell-hook), register now.
# shellcheck disable=SC1091
. "${_enrichviz_comp}"

unset _enrichviz_scripts _enrichviz_root _enrichviz_comp
