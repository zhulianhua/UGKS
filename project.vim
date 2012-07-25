set noautochdir
exec "cd " . escape(expand("<sfile>:p:h"), ' ')

if filereadable("Session.vim")
    source Session.vim
endif

let g:TagCmd='ctags -R --exclude=".git" --exclude=".gitignore" --exclude="*.vim"  --exclude="*.swp"'
let g:SearchPath='**/*.f90'
