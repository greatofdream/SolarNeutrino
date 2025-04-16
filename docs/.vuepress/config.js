import { defaultTheme } from '@vuepress/theme-default'
import { defineUserConfig } from 'vuepress/cli'
import { viteBundler } from '@vuepress/bundler-vite'

export default defineUserConfig({
  lang: 'en-US',
  // change the base url accroding to your setting
  base: '~greatofdream/solar',

  title: 'Solar Neutrino Webpage',
  description: 'Summary and data for the solar neutrino research',

  theme: defaultTheme({
    logo: 'https://ars.els-cdn.com/content/image/1-s2.0-S0146641023000248-gr4.jpg',

    repo: 'https://github.com/greatofdream/SolarNeutrino',

    navbar: ['/', '/SSM', '/SolarNeutrino'],
  }),

  bundler: viteBundler(),
})
