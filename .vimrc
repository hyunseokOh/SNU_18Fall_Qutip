"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" General
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set history=1000
set autoread
set autowrite
set mouse=a

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" VIM user interface
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set ruler
set hlsearch
set incsearch
set showmatch
set number

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Colors and Fonts
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
syntax on
" 컬러스킴
set background=dark
set encoding=utf8

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Text, tab and indent related
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set shiftwidth=4
set tabstop=4
set cindent
set smartindent

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Coding
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set autoindent
set visualbell
set number
set title
set wrap
set linebreak

set nowrapscan
set ignorecase

set backspace=eol,start,indent
set history=100

set nocompatible              " be iMproved, required
filetype off                  " required
	
set rtp+=~/.vim/bundle/Vundle.vim
call vundle#begin()
	Plugin 'VundleVim/Vundle.vim'
	Plugin 'scrooloose/nerdtree'
	Plugin 'Valloric/YouCompleteMe'
	Plugin 'vim-airline/vim-airline'
	Plugin 'pangloss/vim-simplefold'
	Plugin 'nanotech/jellybeans.vim'
	Plugin 'w0rp/ale'
	Plugin 'majutsushi/tagbar'
	Plugin 'scrooloose/nerdcommenter'
call vundle#end()

filetype plugin indent on    " required

"color
colorscheme jellybeans

"NerdTreeOption
let NERDTreeShowHidden=1
let NERDTreeQuitOnOpen=1

"NerdCommenterOption
" Add spaces after comment delimiters by default
let g:NERDSpaceDelims = 1
" Use compact syntax for prettified multi-line comments
let g:NERDCompactSexyComs = 1
" Align line-wise comment delimiters flush left instead of following code
" indentation
let g:NERDDefaultAlign = 'left'
" " Set a language to use its alternate delimiters by default
let g:NERDAltDelims_java = 1
" " Add your own custom formats or override the defaults
let g:NERDCustomDelimiters = { 'c': { 'left': '/**','right': '*/' } }
" " Allow commenting and inverting empty lines (useful when commenting a
" region)
let g:NERDCommentEmptyLines = 1
" " Enable trimming of trailing whitespace when uncommenting
let g:NERDTrimTrailingWhitespace = 1
" " Enable NERDCommenterToggle to check all selected lines is commented or not 
let g:NERDToggleCheckAllLines = 1

"Key mapping
let mapleader=","
map <C-r> :NERDTreeToggle<cr>
map <C-t> :TagbarToggle<cr>
"
	" Brief help
	" :PluginList       - lists configured plugins
	" :PluginInstall    - installs plugins; append `!` to update or just :PluginUpdate
	" :PluginSearch foo - searches for foo; append `!` to refresh local cache
	" :PluginClean      - confirms removal of unused plugins; append `!` to auto-approve removal
	"
	" see :h vundle for more details or wiki for FAQ
	" Put your non-Plugin stuff after this line
