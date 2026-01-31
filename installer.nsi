; ============================
; Schoenflies NSIS Installer
; ============================

!define APP_NAME "ljsim"
!define APP_VERSION "1.1.0"
!define APP_PUBLISHER "Ivo Filot"
!define APP_EXE "ljsim.exe"

OutFile "${APP_NAME}-${APP_VERSION}-Setup.exe"
InstallDir "$PROGRAMFILES64\${APP_NAME}"
InstallDirRegKey HKLM "Software\${APP_NAME}" "InstallDir"

RequestExecutionLevel admin

Page directory
Page instfiles

Section "Install"
    SetOutPath "$INSTDIR"

    ; Copy everything from staging directory
    File /r "dist\${APP_NAME}\*"

    ; Write uninstall information
    WriteRegStr HKLM "Software\${APP_NAME}" "InstallDir" "$INSTDIR"
    WriteUninstaller "$INSTDIR\Uninstall.exe"

    ; Start menu shortcut
    CreateDirectory "$SMPROGRAMS\${APP_NAME}"
    CreateShortcut "$SMPROGRAMS\${APP_NAME}\${APP_NAME}.lnk" "$INSTDIR\${APP_EXE}"

    ; Desktop shortcut
    CreateShortcut "$DESKTOP\${APP_NAME}.lnk" "$INSTDIR\${APP_EXE}"
SectionEnd

Section "Uninstall"
    Delete "$INSTDIR\${APP_EXE}"
    RMDir /r "$INSTDIR"

    Delete "$DESKTOP\${APP_NAME}.lnk"
    RMDir /r "$SMPROGRAMS\${APP_NAME}"

    DeleteRegKey HKLM "Software\${APP_NAME}"
SectionEnd