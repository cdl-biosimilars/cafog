# -*- mode: python -*-

block_cipher = None


a = Analysis(['cafog_gui.py'],
             pathex=[],
             binaries=[],
             datas=[('sample_data', 'sample_data'),
                    ('README.md', '.'),
                    ('LICENSE', '.'),
					('docs/_build', 'docs/_build'),
					('cafog_de.qm', '.')],
             hiddenimports=['pydot'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='cafog_gui',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='cafog_gui')
