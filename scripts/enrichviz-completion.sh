# Tab completion for enrichviz. Sourced from scripts/activate.sh and from
# PROMPT_COMMAND inside `pixi shell` (pixi does not keep `complete` builtins).

if type complete >/dev/null 2>&1 && complete -p enrichviz >/dev/null 2>&1; then
  return 0
fi

_enrichviz() {
  local cur prev track
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  local ora_cmds="gprofiler barplot lollipop dotplot upset ssplot simplify"
  local gsea_cmds="go kegg barplot lollipop dotplot ridgeplot volcano nes nes-compare"

  if [[ ${COMP_CWORD} -eq 1 ]]; then
    COMPREPLY=( $(compgen -W "ora gsea --help" -- "${cur}") )
    return 0
  fi

  if [[ ${COMP_CWORD} -eq 2 ]]; then
    track="${COMP_WORDS[1]}"
    case "${track}" in
      ora) COMPREPLY=( $(compgen -W "${ora_cmds}" -- "${cur}") ) ;;
      gsea) COMPREPLY=( $(compgen -W "${gsea_cmds}" -- "${cur}") ) ;;
    esac
    return 0
  fi

  case "${prev}" in
    --in | --outdir | --rank | --expr)
      COMPREPLY=( $(compgen -f -- "${cur}") )
      ;;
  esac
}

if [ -n "${ZSH_VERSION:-}" ]; then
  autoload -Uz +X bashcompinit 2>/dev/null || autoload -Uz bashcompinit
  bashcompinit 2>/dev/null || true
  if ! typeset -f compdef >/dev/null 2>&1; then
    autoload -Uz compinit
    compinit -C
  fi
fi

if type complete >/dev/null 2>&1; then
  complete -r enrichviz 2>/dev/null || true
  complete -F _enrichviz enrichviz
  complete -F _enrichviz bin/enrichviz
  complete -F _enrichviz ./bin/enrichviz
fi
