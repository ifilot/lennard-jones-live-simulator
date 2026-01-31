; ============================
; Schoenflies NSIS Installer
; ============================

!define APP_NAME "ljsim"
!define APP_VERSION "1.1.0"
!define APP_PUBLISHER "Ivo Filot"
!define APP_EXE "ljsim.exe"

; --- Installer / Uninstaller icon ---
Icon "assets\icon\ljsim.ico"

OutFile "${APP_NAME}-${APP_VERSION}-setup.exe"
InstallDir "$PROGRAMFILES64\${APP_NAME}"
InstallDirRegKey HKLM "Software\${APP_NAME}" "InstallDir"

RequestExecutionLevel admin

Page directory
Page instfiles

Section "Install"
    SetOutPath "$INSTDIR"

    ; Copy shortcut icon into installation directory
    File "assets\icon\ljsim.ico"

    ; Copy everything from staging directory
    File /r "dist\${APP_NAME}\*"

    ; Write uninstall information
    WriteRegStr HKLM "Software\${APP_NAME}" "InstallDir" "$INSTDIR"
    WriteUninstaller "$INSTDIR\Uninstall.exe"

    ; Start menu shortcut with custom icon
    CreateDirectory "$SMPROGRAMS\${APP_NAME}"
    CreateShortcut "$SMPROGRAMS\${APP_NAME}\${APP_NAME}.lnk" "$INSTDIR\${APP_EXE}" "" "$INSTDIR\ljsim.ico"

    ; Desktop shortcut with custom icon
    CreateShortcut "$DESKTOP\${APP_NAME}.lnk" "$INSTDIR\${APP_EXE}" "" "$INSTDIR\ljsim.ico"
SectionEnd

Section "Uninstall"
    Delete "$INSTDIR\${APP_EXE}"
    Delete "$INSTDIR\ljsim.ico"
    Delete "$INSTDIR\Uninstall.exe"
    RMDir /r "$INSTDIR"

    Delete "$DESKTOP\${APP_NAME}.lnk"
    RMDir /r "$SMPROGRAMS\${APP_NAME}"

    DeleteRegKey HKLM "Software\${APP_NAME}"
SectionEnd
